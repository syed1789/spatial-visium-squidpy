FROM python:3.13-slim
# g++ is needed to build mahotas (a Squidpy dependency) where no wheel exists, e.g. linux/arm64
RUN apt-get update && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY analyze_visium.py .
CMD ["python", "analyze_visium.py", "--outdir", "results"]
