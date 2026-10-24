FROM --platform=linux/amd64 python:3.12.15-slim@sha256:05cda9777409a9c3ffddd94a4c476b79f0769a0b4857f0c7ed9226b6800b0d6f
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
WORKDIR /app
COPY requirements.lock.txt /tmp/runtime-requirements.lock.txt
RUN python -m pip install --no-cache-dir --only-binary=:all: --no-deps --require-hashes -r /tmp/runtime-requirements.lock.txt && python -m pip check
RUN python -c "import numpy as np; from plasmapy.dispersion import plasma_dispersion_func as Z; assert abs(Z(0)-1j*np.sqrt(np.pi))<1e-14"
COPY src/ /app/src/
COPY runner/ /runner/
ENV PYTHONPATH=/app/src
