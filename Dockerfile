FROM python:3.13.14-trixie AS development

COPY --from=ghcr.io/astral-sh/uv:0.11.23 /uv /uvx /usr/local/bin/

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPYCACHEPREFIX=/tmp/pycache \
    UV_LINK_MODE=copy

WORKDIR /workspaces/vcf_annotator

FROM development AS testing

COPY pyproject.toml uv.lock ./
RUN uv sync --extra test --frozen

COPY . .
RUN uv run pytest

FROM python:3.13.14-slim-trixie AS package

COPY --from=ghcr.io/astral-sh/uv:0.11.23 /uv /uvx /usr/local/bin/

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPYCACHEPREFIX=/tmp/pycache \
    UV_LINK_MODE=copy

WORKDIR /app

COPY pyproject.toml uv.lock ./
COPY vcf_annotator ./vcf_annotator

RUN uv build --wheel --out-dir /dist

FROM python:3.13.14-slim-trixie AS build

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPYCACHEPREFIX=/tmp/pycache \
    PATH=/opt/venv/bin:$PATH

WORKDIR /app

COPY --from=package /dist/*.whl /tmp/
RUN python -m venv /opt/venv \
    && pip install --no-cache-dir /tmp/*.whl \
    && rm -rf /tmp/*.whl /root/.cache

ENTRYPOINT ["python", "-m", "vcf_annotator.cli"]
