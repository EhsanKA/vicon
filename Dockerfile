FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Copy environment and package files
COPY environment.yaml .
COPY pyproject.toml .
COPY setup.py .

# Create conda environment
RUN conda env create -f environment.yaml

# Install the package in development mode
RUN conda run -n vicon pip install -e .

# Copy the rest of the package
COPY . .

# Add wrapper script directory to PATH
ENV PATH="/app/scripts:${PATH}"

# Default command
CMD ["bash"]
