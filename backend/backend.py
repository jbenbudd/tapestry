from typing import Optional
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from backend.pubmed import search_and_fetch

app = FastAPI(
    title="Tapestry Backend",
    description="Backend for Tapestry: A retrieval-augmented generation (RAG) engine for academic research",
    version="0.1.0"
)

# Configure CORS so that your React front end (or any other client) can access the API.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # For development, '*' is fine; later, restrict to your front end's origin.
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/api/search")
async def search_articles(query: str, max_results: int = 10, start_date: Optional[str] = None, end_date: Optional[str] = None):
    """
    Search for PubMed articles based on the query.
    Returns a list of articles with id, title, and abstract.
    """
    try:
        articles = search_and_fetch(query, max_results, start_date, end_date)
        if not articles:
            raise HTTPException(status_code=404, detail="No articles found.")
        return articles
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == '__main__':
    import uvicorn
    uvicorn.run("backend.backend:app", host="0.0.0.0", port=8000, reload=True)