FROM python:3.12-slim-bookworm

# Keep Python fast/noisy in containers
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# System deps for GDAL/GEOS/PROJ + build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ make \
    gdal-bin libgdal-dev \
    libproj-dev proj-bin \
    libgeos-dev \
    libspatialindex-dev \
    libtiff-dev libjpeg-dev libpng-dev \
    libcurl4-openssl-dev \
    procps \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Help pip find GDAL headers (Fiona/Rasterio look for these)
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal \
    C_INCLUDE_PATH=/usr/include/gdal \
    GDAL_CONFIG=/usr/bin/gdal-config

# Python deps
RUN pip install --no-cache-dir \
    numpy pandas scipy matplotlib plotly \
    shapely \
    pyproj \
    fiona \
    rasterio \
    geopandas \
    folium \
    pykrige

WORKDIR /app
