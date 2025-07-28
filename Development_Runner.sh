#!/bin/bash
# filepath: /fast/AG_Ohler/ekarimi/projects/vicon/run_web.sh

echo "🔧 Starting VICON Web Interface in development mode..."

# Create directories
mkdir -p web_app/{jobs,static/{css,js},templates}

# Install dependencies
if [ ! -f "venv/bin/activate" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
fi

source venv/bin/activate

echo "📥 Installing dependencies..."
pip install -r requirements.txt
pip install -r requirements-web.txt

# Set environment variables
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
export ENVIRONMENT="development"

# Run the application
echo "🚀 Starting web server..."
echo "🌐 Access the application at: http://localhost:8000"
echo "📚 API documentation at: http://localhost:8000/api/docs"

uvicorn web_app.main:app --host 0.0.0.0 --port 8000 --reload