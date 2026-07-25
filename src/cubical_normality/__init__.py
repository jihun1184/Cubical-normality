"""Embedded 3D cubical-normality checker used by the accompanying paper."""

from .checker import (
    CERTIFICATE_BACKEND,
    CERTIFICATE_SCHEMA_VERSION,
    InputValidationError,
    build_certificate,
    build_certificate_cli,
    check_embedded_cubical_normality,
)

__all__ = [
    "CERTIFICATE_BACKEND",
    "CERTIFICATE_SCHEMA_VERSION",
    "InputValidationError",
    "build_certificate",
    "build_certificate_cli",
    "check_embedded_cubical_normality",
]
