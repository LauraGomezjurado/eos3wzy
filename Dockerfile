FROM bentoml/model-server:0.11.0-py39
MAINTAINER ersilia

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    gcc \
    g++ \
    make \
    libopenblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --upgrade pip
RUN pip install pandas numpy scipy
RUN pip install rdkit
RUN pip install pytorch-lightning
RUN pip install torch-geometric
RUN pip install openbabel
RUN pip install multiprocessing-logging

# Install QupKake
RUN pip install git+https://github.com/Shualdon/QupKake.git

# Install xtb (GFN2-xTB)
RUN wget https://github.com/grimme-lab/xtb/releases/download/v6.4.1/xtb-6.4.1-linux-x86_64.tar.xz && \
    tar -xf xtb-6.4.1-linux-x86_64.tar.xz && \
    mv xtb-6.4.1 /opt/xtb && \
    rm xtb-6.4.1-linux-x86_64.tar.xz

# Add xtb to PATH
ENV PATH="/opt/xtb/bin:${PATH}"
ENV XTBPATH="/opt/xtb"

# Create cache directory
RUN mkdir -p /app/cache

WORKDIR /repo
COPY . /repo

# Make CLI executable
RUN chmod +x /repo/src/cli.py

# Set default command
CMD ["python", "/repo/src/cli.py", "info"]
