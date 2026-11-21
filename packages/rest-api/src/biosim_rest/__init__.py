# SAVE TO: packages\rest-api\src\biosim_rest\__init__.py
"""
biosim_rest -- HTTP REST layer over biosim_api.

Usage:
    from biosim_rest import create_app
    app = create_app()

    # Run with: uvicorn biosim_rest:create_app --factory --reload
"""

from .app import create_app

__all__ = ["create_app"]
__version__ = "0.1.0"
