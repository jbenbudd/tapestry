# backend/arxiv.py

import feedparser


def search_arxiv(query, max_results=10):
    """Search arXiv for a given query and return feed entries."""
    base_url = "http://export.arxiv.org/api/query?"
    # arXiv API accepts a search query in the format 'search_query=all:QUERY'
    search_query = f"search_query=all:{query}&start=0&max_results={max_results}"
    url = base_url + search_query
    feed = feedparser.parse(url)
    return feed.entries


def fetch_arxiv_details(query, max_results=10):
    """Fetch and parse arXiv article details based on a query."""
    entries = search_arxiv(query, max_results)
    articles = []
    for entry in entries:
        title = entry.get('title', 'No title available').strip().replace('\n', ' ')
        abstract = entry.get('summary', 'No abstract available').strip().replace('\n', ' ')
        authors = [author.name for author in entry.get('authors', [])]
        authors_str = ", ".join(authors) if authors else "No authors available"
        published = entry.get('published', 'No publication date available')
        link = entry.get('id', 'No link available')

        articles.append({
            "id": entry.get('id'),
            "title": title,
            "abstract": abstract,
            "authors": authors_str,
            "journal": "arXiv",
            "pub_date": published,
            "link": link
        })
    return articles


def search_and_fetch(query, max_results=10):
    """Search and fetch articles from arXiv based on a query."""
    return fetch_arxiv_details(query, max_results) 