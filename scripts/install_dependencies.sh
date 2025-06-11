#!/bin/bash
# This script installs the required tools: Docker, Java, Nextflow and downloads the SILVA reference database.

set -e

echo "Updating system and installing dependencies..."

# Update package lists
sudo apt update

# Install Docker
echo "Installing Docker..."
sudo apt install docker.io docker-compose

# Add current user to the Docker group
sudo usermod -aG docker $USER

# Install Java (required by Nextflow)
echo "Installing Java..."
sudo apt install default-jre

# Install Nextflow
echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin

# Download SILVA reference database (v138.1, seed release)
echo "Downloading SILVA v138.1 (seed release)..."
mkdir -p ~/databases
cd ~/databases
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz
tar xfv silva.seed_v138_1.tgz
rm silva.seed_v138_1.tgz

echo "All dependencies installed successfully."
