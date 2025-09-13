import uvicorn
from fastapi import FastAPI, Request, HTTPException
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, field_validator
from tasks.app import run_tool
from enum import Enum
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Sequence Alignment API",
    description="API for running sequence alignment tools (MAFFT, VSEARCH, UCLUST)",
    version="1.0.0",
)

class AlignTool(str, Enum):
    MAFFT = "mafft"
    VSEARCH = "vsearch"
    UCLUST = "uclust"


class RunRequest(BaseModel):
    dir_name: str = Field(..., description="directory name")
    base_name: str = Field(..., description="input file name")
    align_tool: str = Field(
        ..., description="Alignment tool to use (mafft, vsearch, or uclust)"
    )
    options: str = Field(..., description="Command line options for the tool")


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    logger.error(f"Validation error for {request.method} {request.url}: {exc.errors()}")
    return JSONResponse(
        status_code=422,
        content={
            "error": "Request validation failed",
            "details": exc.errors(),
            "expected_format": {
                "align_tool": "One of: mafft, vsearch, uclust",
                "options": "Command line options as string",
            },
        },
    )


@app.post("/align")
async def trigger_run(request: RunRequest):
    logger.info(
        f"Received alignment request for file: {request.base_name}, tool: {request.align_tool}, options: {request.options}"
    )

    try:
        task = run_tool.delay(
            dir_name=request.dir_name,
            base_name=request.base_name,
            tool=request.align_tool,
            options=request.options,
        )
        logger.info(f"Task queued with ID: {task.id}")
        return {"status": "PENDING", "task_id": task.id}
    except Exception as e:
        logger.error(f"Error queuing task: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to queue alignment task: {str(e)}"
        )


@app.get("/")
async def root():
    return {
        "message": "Sequence Alignment API",
        "endpoints": {
            "POST /align": "Trigger alignment job",
            "GET /docs": "API documentation",
            "GET /health": "Health check",
        },
    }


@app.get("/health")
async def health_check():
    return {"status": "healthy", "service": "sequence-alignment-api"}

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)