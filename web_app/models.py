from pydantic import BaseModel, EmailStr, Field, validator
from typing import Optional, Dict, Any
from datetime import datetime
from enum import Enum

class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class PipelineConfig(BaseModel):
    """Pipeline configuration model"""
    project_path: str = ""
    virus_name: str = Field(..., min_length=1, max_length=100)
    input_sample: str = "data/sample.fasta"
    input_reference: str = "data/reference.fasta"
    email: EmailStr
    kmer_size: int = Field(150, ge=1, le=1000)
    threshold: int = Field(147, ge=1)
    only_valid_kmers: bool = True
    l_gene_start: int = -1
    l_gene_end: int = Field(16000, ge=0)
    coverage_ratio: float = Field(0.6, ge=0.0, le=1.0)
    sort_by_mismatches: bool = True
    min_year: int = Field(2020, ge=1900, le=2030)
    threshold_ratio: float = Field(0.01, ge=0.0, le=1.0)
    drop_old_samples: bool = False
    drop_mischar_samples: bool = True

    @validator('threshold')
    def validate_threshold(cls, v, values):
        if 'kmer_size' in values and v >= values['kmer_size']:
            raise ValueError('Threshold must be less than kmer_size')
        return v

class JobInfo(BaseModel):
    """Job information model"""
    job_id: str
    status: JobStatus
    created_at: datetime
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    config: Dict[str, Any]
    results_available: bool = False

class JobResponse(BaseModel):
    """Job submission response"""
    job_id: str
    status: str
    message: str

class JobStatusResponse(BaseModel):
    """Job status response"""
    job_id: str
    status: JobStatus
    created_at: datetime
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    results_available: bool = False
    progress: Optional[str] = None