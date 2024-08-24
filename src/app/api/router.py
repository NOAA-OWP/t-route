from fastapi import APIRouter

from app.api.routes import v4_routing

api_router = APIRouter()
api_router.include_router(v4_routing.router, prefix="/flow_routing/v4", tags=["v4 Flow Routing"])
