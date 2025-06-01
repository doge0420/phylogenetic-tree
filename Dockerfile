# ─── Stage 1: Base image ─────────────────────────────────────
FROM python:3.10-slim

# ─── Stage 2: Set working directory ──────────────────────────
WORKDIR /app

# ─── Stage 3: Install Python dependencies ────────────────────
RUN pip install --no-cache-dir numpy matplotlib

# ─── Stage 4: Prepare output folder and declare as volume ────
RUN mkdir -p /app/output
VOLUME /app/output

# ─── Stage 5: Copy your updated upgma.py ──────────────────────
COPY upgma.py /app/upgma.py

# ─── Stage 6: Default command ─────────────────────────────────
# By default, run in “sequence mode” (no --distances). 
# If you want matrix mode, you’ll override this in `podman run`.
CMD ["python", "upgma.py"]
