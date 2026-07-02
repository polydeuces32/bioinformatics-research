# Build context is the repository root (see cloudflare/wrangler.jsonc
# image_build_context). This image only serves the dashboard — it does not
# need pydeseq2/gseapy/lifelines/GEOparse, so it uses its own lean
# requirements.txt rather than environment.yml.
FROM python:3.11-slim

WORKDIR /app

COPY dashboard/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY dashboard/app.py .

# Bakes in whatever results/tables/ contains at build time (a point-in-time
# snapshot, since results/ is gitignored and generated locally). Rebuild and
# redeploy the container to publish an updated snapshot — this dashboard
# does not run the analysis pipeline itself.
COPY results/tables ./results/tables

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.headless=true"]
