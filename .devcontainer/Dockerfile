# Start from the Python 3.11 slim image as a base
FROM python:3.11-slim

# Install prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates gnupg \
    && rm -rf /var/lib/apt/lists/*

# Install Java (OpenJDK). Nextflow supports Java 8+, commonly Java 11 or 17.
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-11-jdk \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME environment variable for convenience
ENV JAVA_HOME="/usr/lib/jvm/java-11-openjdk-amd64"
ENV PATH="$JAVA_HOME/bin:$PATH"

# Install Python packages
RUN pip install --no-cache-dir pandas numpy

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

# Create a non-root user (named 'vscode') to match devcontainer defaults
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

USER $USERNAME
WORKDIR /workspaces
