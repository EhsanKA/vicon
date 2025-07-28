#!/bin/bash
# filepath: /fast/AG_Ohler/ekarimi/projects/vicon/deploy.sh

set -e

echo "🚀 Deploying VICON Web Interface..."

# Configuration
DOMAIN="your-domain.com"
EMAIL="your-email@domain.com"

# Create necessary directories
mkdir -p web_app/{jobs,static/{css,js},templates}
mkdir -p ssl

# Generate SSL certificates (using Let's Encrypt)
if [ ! -f "ssl/cert.pem" ]; then
    echo "📜 Generating SSL certificates..."
    # Install certbot if not present
    if ! command -v certbot &> /dev/null; then
        sudo apt update
        sudo apt install -y certbot
    fi
    
    # Get certificates
    sudo certbot certonly --standalone -d $DOMAIN --email $EMAIL --agree-tos --non-interactive
    
    # Copy certificates
    sudo cp /etc/letsencrypt/live/$DOMAIN/fullchain.pem ssl/cert.pem
    sudo cp /etc/letsencrypt/live/$DOMAIN/privkey.pem ssl/key.pem
    sudo chown $USER:$USER ssl/*.pem
fi

# Build and deploy with Docker Compose
echo "🐳 Building and starting containers..."
docker-compose down --remove-orphans
docker-compose build --no-cache
docker-compose up -d

# Wait for services to be ready
echo "⏳ Waiting for services to start..."
sleep 30

# Health check
echo "🏥 Performing health check..."
if curl -f http://localhost/health > /dev/null 2>&1; then
    echo "✅ VICON Web Interface is running!"
    echo "🌐 Access it at: https://$DOMAIN"
else
    echo "❌ Health check failed!"
    docker-compose logs
    exit 1
fi

# Setup log rotation
echo "📝 Setting up log rotation..."
sudo tee /etc/logrotate.d/vicon << EOF
/var/lib/docker/containers/*/*.log {
    rotate 7
    daily
    compress
    delaycompress
    missingok
    notifempty
    create 0644 root root
}
EOF

echo "🎉 Deployment complete!"
echo "📊 Monitor logs with: docker-compose logs -f"
echo "🔄 Update with: git pull && docker-compose up -d --build"