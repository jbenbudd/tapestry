import React, { useState } from 'react';

interface Article {
  id: string;
  title: string;
  abstract: string;
}

const App: React.FC = () => {
  // State for the search query, results, selected articles, and the currently previewed article.
  const [query, setQuery] = useState<string>('');
  const [results, setResults] = useState<Article[]>([]);
  const [selectedArticles, setSelectedArticles] = useState<Article[]>([]);
  const [previewArticleId, setPreviewArticleId] = useState<string | null>(null);

  // Handler for searching articles
  const handleSearch = async () => {
    if (!query) return;
    try {
      // Replace the URL with your actual backend endpoint
      const response = await fetch(`http://127.0.0.1:8000/api/search?query=${encodeURIComponent(query)}`);
      if (!response.ok) {
        console.error('Search failed.');
        return;
      }
      const data: Article[] = await response.json();
      setResults(data);
      setPreviewArticleId(null);
    } catch (error) {
      console.error(error);
    }
  };

  // Handler to add an article to the selected list
  const handleAdd = (article: Article) => {
    if (!selectedArticles.find(a => a.id === article.id)) {
      setSelectedArticles([...selectedArticles, article]);
    }
  };

  // Handler to remove an article from the selected list
  const handleRemove = (articleId: string) => {
    setSelectedArticles(selectedArticles.filter(article => article.id !== articleId));
  };

  return (
    <div className="min-h-screen bg-gray-100 p-4">
      <h1 className="text-3xl font-bold text-center mb-6">
        Tapestry: A retrieval-augmented generation engine for academic research
      </h1>
      <div className="max-w-4xl mx-auto">
        {/* Search Section */}
        <div className="mb-6">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Search PubMed articles..."
            className="w-full p-2 border border-gray-300 rounded-md"
          />
          <button
            onClick={handleSearch}
            className="mt-2 px-4 py-2 bg-blue-500 text-white rounded-md hover:bg-blue-600"
          >
            Search
          </button>
        </div>
        <div className="grid grid-cols-2 gap-4">
          {/* Search Results Section */}
          <div className="bg-white p-4 rounded-md shadow-md">
            <h2 className="text-xl font-semibold mb-4">Search Results</h2>
            {results.length === 0 ? (
              <p className="text-gray-500">No results yet.</p>
            ) : (
              <ul>
                {results.map(article => (
                  <li key={article.id} className="border-b border-gray-200 py-3">
                    <div className="flex justify-between items-center">
                      <div>
                        <p className="font-bold">{article.title}</p>
                        {/* Only show abstract when this article is set to be previewed */}
                        {previewArticleId === article.id && (
                          <p className="text-sm text-gray-700 mt-1">{article.abstract}</p>
                        )}
                      </div>
                      <div className="flex flex-col items-end space-y-1">
                        <button
                          onClick={() => setPreviewArticleId(
                            previewArticleId === article.id ? null : article.id
                          )}
                          className="text-sm text-blue-500 hover:underline"
                        >
                          {previewArticleId === article.id ? 'Hide' : 'Preview'}
                        </button>
                        <button
                          onClick={() => handleAdd(article)}
                          className="text-sm text-green-500 hover:underline"
                        >
                          Add
                        </button>
                      </div>
                    </div>
                  </li>
                ))}
              </ul>
            )}
          </div>
          {/* Selected Articles Section */}
          <div className="bg-white p-4 rounded-md shadow-md">
            <h2 className="text-xl font-semibold mb-4">Selected Articles</h2>
            {selectedArticles.length === 0 ? (
              <p className="text-gray-500">No articles selected.</p>
            ) : (
              <ul>
                {selectedArticles.map(article => (
                  <li key={article.id} className="border-b border-gray-200 py-3 flex justify-between items-center">
                    <p className="font-bold">{article.title}</p>
                    <button
                      onClick={() => handleRemove(article.id)}
                      className="text-sm text-red-500 hover:underline"
                    >
                      Remove
                    </button>
                  </li>
                ))}
              </ul>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default App;