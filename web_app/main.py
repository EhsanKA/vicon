from fastapi import FastAPI, File, UploadFile, Form, BackgroundTasks, Request, HTTPException, Depends
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import HTTPBasic, HTTPBasicCredentials
import yaml
import os
import uuid
import asyncio
import subprocess
import shutil
import logging
from typing import Optional, Dict, List
import json
from datetime import datetime
from pathlib import Path

from .models import PipelineConfig, JobInfo, JobResponse, JobStatusResponse, JobStatus

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="VICON Pipeline Web Interface", 
    description="Viral Consensus Pipeline - Web-based interface for running VICON analysis",
    version="1.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc"
)

# CORS middleware for cross-origin requests
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Basic authentication (optional)
security = HTTPBasic()

# Ensure directories exist
WEB_APP_DIR = Path(__file__).parent
JOBS_DIR = WEB_APP_DIR / "jobs"
STATIC_DIR = WEB_APP_DIR / "static"
TEMPLATES_DIR = WEB_APP_DIR / "templates"

JOBS_DIR.mkdir(exist_ok=True)
STATIC_DIR.mkdir(exist_ok=True)
TEMPLATES_DIR.mkdir(exist_ok=True)

# Mount static files and templates
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))

# In-memory job storage (use Redis/Database for production)
jobs: Dict[str, JobInfo] = {}

# Configuration limits
MAX_FILE_SIZE = 500 * 1024 * 1024  # 500MB
ALLOWED_EXTENSIONS = {'.fasta', '.fa', '.fas'}

def validate_file(file: UploadFile) -> None:
    """Validate uploaded file"""
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file provided")
    
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid file type. Allowed: {', '.join(ALLOWED_EXTENSIONS)}"
        )

def get_current_user(credentials: HTTPBasicCredentials = Depends(security)):
    """Basic authentication (optional)"""
    # Implement your authentication logic here
    return credentials.username

@app.get("/", response_class=HTMLResponse)
async def main_form(request: Request):
    """Main pipeline configuration form"""
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "timestamp": datetime.now()}

@app.post("/api/run-pipeline/", response_model=JobResponse)
async def run_pipeline(
    background_tasks: BackgroundTasks,
    virus_name: str = Form(...),
    email: str = Form(...),
    sample_file: UploadFile = File(...),
    reference_file: UploadFile = File(...),
    kmer_size: int = Form(150),
    threshold: int = Form(147),
    only_valid_kmers: bool = Form(True),
    l_gene_start: int = Form(-1),
    l_gene_end: int = Form(16000),
    coverage_ratio: float = Form(0.6),
    sort_by_mismatches: bool = Form(True),
    min_year: int = Form(2020),
    threshold_ratio: float = Form(0.01),
    drop_old_samples: bool = Form(False),
    drop_mischar_samples: bool = Form(True)
):
    """Submit a new pipeline job"""
    
    # Validate files
    validate_file(sample_file)
    validate_file(reference_file)
    
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    
    data_dir = job_dir / "data"
    data_dir.mkdir(exist_ok=True)
    
    try:
        # Save uploaded files
        sample_path = data_dir / "sample.fasta"
        reference_path = data_dir / "reference.fasta"
        
        with open(sample_path, "wb") as f:
            content = await sample_file.read()
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(status_code=413, detail="File too large")
            f.write(content)
            
        with open(reference_path, "wb") as f:
            content = await reference_file.read()
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(status_code=413, detail="File too large")
            f.write(content)
        
        # Create configuration
        config_data = {
            "project_path": str(job_dir.absolute()),
            "virus_name": virus_name,
            "input_sample": "data/sample.fasta",
            "input_reference": "data/reference.fasta",
            "email": email,
            "kmer_size": kmer_size,
            "threshold": threshold,
            "only_valid_kmers": only_valid_kmers,
            "l_gene_start": l_gene_start,
            "l_gene_end": l_gene_end,
            "coverage_ratio": coverage_ratio,
            "sort_by_mismatches": sort_by_mismatches,
            "min_year": min_year,
            "threshold_ratio": threshold_ratio,
            "drop_old_samples": drop_old_samples,
            "drop_mischar_samples": drop_mischar_samples
        }
        
        # Validate configuration
        config = PipelineConfig(**config_data)
        
        # Save config file
        config_path = job_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_data, f, default_flow_style=False)
        
        # Store job info
        jobs[job_id] = JobInfo(
            job_id=job_id,
            status=JobStatus.PENDING,
            created_at=datetime.now(),
            config=config_data
        )
        
        # Start background task
        background_tasks.add_task(run_vicon_pipeline, job_id, str(config_path))
        
        logger.info(f"Job {job_id} submitted for virus: {virus_name}")
        
        return JobResponse(
            job_id=job_id,
            status="pending",
            message="Pipeline job submitted successfully"
        )
        
    except Exception as e:
        # Cleanup on error
        if job_dir.exists():
            shutil.rmtree(job_dir)
        if job_id in jobs:
            del jobs[job_id]
        
        logger.error(f"Error submitting job: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to submit job: {str(e)}")

async def run_vicon_pipeline(job_id: str, config_path: str):
    """Run the VICON pipeline in background"""
    try:
        jobs[job_id].status = JobStatus.RUNNING
        logger.info(f"Starting pipeline for job {job_id}")
        
        # Get project root directory
        project_root = Path(__file__).parent.parent
        
        # Prepare command
        cmd = [
            "python", "-m", "vicon.run_pipeline", 
            "--config", config_path
        ]
        
        # Run the pipeline
        process = await asyncio.create_subprocess_exec(
            *cmd,
            cwd=str(project_root),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        stdout, stderr = await process.communicate()
        
        if process.returncode == 0:
            jobs[job_id].status = JobStatus.COMPLETED
            jobs[job_id].completed_at = datetime.now()
            
            # Check if results are available
            job_dir = JOBS_DIR / job_id
            results_dir = job_dir / "results" / jobs[job_id].config["virus_name"]
            jobs[job_id].results_available = results_dir.exists()
            
            logger.info(f"Pipeline completed successfully for job {job_id}")
        else:
            jobs[job_id].status = JobStatus.FAILED
            error_msg = stderr.decode() if stderr else "Unknown error occurred"
            jobs[job_id].error_message = error_msg
            jobs[job_id].completed_at = datetime.now()
            
            logger.error(f"Pipeline failed for job {job_id}: {error_msg}")
            
    except Exception as e:
        jobs[job_id].status = JobStatus.FAILED
        jobs[job_id].error_message = str(e)
        jobs[job_id].completed_at = datetime.now()
        
        logger.error(f"Exception in pipeline for job {job_id}: {str(e)}")

@app.get("/api/status/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """Get job status"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Check for results if completed
    results_available = False
    if job.status == JobStatus.COMPLETED:
        job_dir = JOBS_DIR / job_id
        results_dir = job_dir / "results" / job.config["virus_name"]
        results_available = results_dir.exists()
    
    return JobStatusResponse(
        job_id=job_id,
        status=job.status,
        created_at=job.created_at,
        completed_at=job.completed_at,
        error_message=job.error_message,
        results_available=results_available
    )

@app.get("/api/download/{job_id}/{file_type}")
async def download_results(job_id: str, file_type: str):
    """Download result files"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    job_dir = JOBS_DIR / job_id
    results_dir = job_dir / "results" / job.config["virus_name"]
    
    file_paths = {
        "kmer1": results_dir / "kmer1.csv",
        "kmer2": results_dir / "kmer2.csv",
        "logs": results_dir / "logs" / "pipeline.log",
        "config": job_dir / "config.yaml",
        "plot": results_dir / f"{job.config['virus_name']}_coverage_plot.png"
    }
    
    if file_type not in file_paths:
        raise HTTPException(status_code=400, detail="Invalid file type")
    
    file_path = file_paths[file_type]
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    return FileResponse(
        str(file_path),
        filename=f"{job.config['virus_name']}_{file_type}.{file_path.suffix}",
        media_type='application/octet-stream'
    )

@app.get("/api/jobs/", response_model=List[JobStatusResponse])
async def list_jobs():
    """List all jobs"""
    return [
        JobStatusResponse(
            job_id=job.job_id,
            status=job.status,
            created_at=job.created_at,
            completed_at=job.completed_at,
            error_message=job.error_message,
            results_available=job.results_available
        ) for job in jobs.values()
    ]

@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job and its files"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # Remove job directory
    job_dir = JOBS_DIR / job_id
    if job_dir.exists():
        shutil.rmtree(job_dir)
    
    # Remove from memory
    del jobs[job_id]
    
    logger.info(f"Job {job_id} deleted")
    return {"message": "Job deleted successfully"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)