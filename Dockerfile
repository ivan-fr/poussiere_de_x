FROM ubuntu:22.04

# Avoid tzdata interactive prompt during apt installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for Lean, elan, and mathlib
RUN apt-get update && apt-get install -y \
    curl \
    git \
    make \
    cmake \
    python3 \
    python3-pip \
    tar \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Add a non-root user (standard for dev containers)
ARG USERNAME=lean_user
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Switch to the new user
USER $USERNAME
ENV HOME=/home/$USERNAME

# Install elan (the official Lean version manager, like rustup for rust)
ENV ELAN_HOME=$HOME/.elan
ENV PATH=$ELAN_HOME/bin:$PATH

RUN curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none

# Set the working directory
WORKDIR /workspace

# Keep container running continuously (useful for docker-compose and dev containers)
CMD ["sleep", "infinity"]
