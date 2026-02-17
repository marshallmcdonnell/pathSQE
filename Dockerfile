# Use a Python-based image with conda pre-installed
FROM continuumio/miniconda3:latest

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Set working directory
WORKDIR /app

# Install pixi using conda
RUN conda install -y -c conda-forge pixi && \
    conda clean --all -y

# Copy project files
COPY . .

# Install dependencies and package using pixi
RUN pixi install

# Set entrypoint to use pixi's Python environment
ENV PATH="/app/.pixi/envs/default/bin:$PATH"

# Default command
CMD ["/bin/bash"]
