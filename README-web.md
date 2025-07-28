# VICON Web Interface

A web-based interface for the VICON viral consensus pipeline.

## Features

- 🌐 Web-based form for pipeline configuration
- 📁 File upload for FASTA samples and references  
- ⚡ Real-time job status monitoring
- 📊 Result visualization and download
- 🐳 Docker containerization
- 🔒 SSL/HTTPS support
- 🎛️ REST API for programmatic access

## Quick Start

### To Make Scripts Executable
```bash
chmod +x run_web.sh deploy.sh
```

### Development Mode
```bash
./run_web.sh
```

### Production Deployment
```bash
# Configure your domain in deploy.sh
./deploy.sh
```

### Using Docker Compose
```bash
docker-compose up -d
```

## API Endpoints

- `GET /` - Main web interface
- `POST /api/run-pipeline/` - Submit pipeline job
- `GET /api/status/{job_id}` - Get job status
- `GET /api/download/{job_id}/{file_type}` - Download results
- `GET /api/jobs/` - List all jobs
- `DELETE /api/jobs/{job_id}` - Delete job

## Configuration

The web interface accepts the same parameters as the CLI version:

- **virus_name**: Identifier for your analysis
- **email**: Contact email (required for some tools)
- **sample_file**: FASTA file with sample sequences
- **reference_file**: FASTA file with reference sequence
- **kmer_size**: Size of k-mers (default: 150)
- **threshold**: Matching threshold (default: 147)
- **l_gene_start/end**: Gene region boundaries
- **coverage_ratio**: Minimum coverage ratio
- **Additional parameters**: See the web form for all options

## File Formats

- **Input**: FASTA format (.fasta, .fa)
- **Output**: CSV files, logs, plots

## Security

- File upload validation
- Rate limiting
- SSL/HTTPS encryption
- Input sanitization
- Container isolation

## Monitoring

- Health check endpoint: `/health`
- Real-time job status updates
- Comprehensive logging
- Error handling and reporting

## Troubleshooting

### Common Issues

1. **File upload fails**: Check file size (<500MB) and format (.fasta/.fa)
2. **Job stuck in pending**: Check Docker containers and logs
3. **SSL certificate issues**: Ensure domain points to server, check certbot

### Logs

```bash
# Application logs
docker-compose logs vicon-web

# Nginx logs  
docker-compose logs nginx

# All services
docker-compose logs -f
```

## Support

For issues and questions, check the logs and ensure all dependencies are installed correctly.