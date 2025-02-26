import React, { useState } from 'react';
import DatePicker from "react-datepicker";
import "react-datepicker/dist/react-datepicker.css";

interface Article {
  id: string;
  title: string;
  abstract: string;
  authors?: string;
  journal?: string;
  pub_date?: string;
  link?: string;
}

const App: React.FC = () => {
  const [query, setQuery] = useState<string>('');
  const [startDate, setStartDate] = useState<Date | null>(null);
  const [endDate, setEndDate] = useState<Date | null>(null);
  const [results, setResults] = useState<Article[]>([]);
  const [selectedArticles, setSelectedArticles] = useState<Article[]>([]);
  const [previewArticleId, setPreviewArticleId] = useState<string | null>(null);
  const [isFilterOpen, setIsFilterOpen] = useState<boolean>(false);

  const handleSearch = async () => {
    if (!query) return;
    try {
      // Format dates as YYYY-MM-DD if they exist
      const formattedStartDate = startDate ? startDate.toISOString().split('T')[0] : "";
      const formattedEndDate = endDate ? endDate.toISOString().split('T')[0] : "";
      const response = await fetch(
        `http://127.0.0.1:8000/api/search?query=${encodeURIComponent(query)}&start_date=${encodeURIComponent(formattedStartDate)}&end_date=${encodeURIComponent(formattedEndDate)}`
      );
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

  const handleAdd = (article: Article) => {
    if (!selectedArticles.find(a => a.id === article.id)) {
      setSelectedArticles([...selectedArticles, article]);
    }
  };

  const handleRemove = (articleId: string) => {
    setSelectedArticles(selectedArticles.filter(article => article.id !== articleId));
  };

  return (
    <div className="min-h-screen bg-gray-100 p-4">
      <h1 className="text-3xl font-bold text-center mb-6">
        Tapestry: A retrieval-augmented generation engine for academic research
      </h1>
      <div className="max-w-4xl mx-auto">
        {/* Main Search Bar (always visible) */}
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

        {/* Filter Menu Toggle */}
        <button
          onClick={() => setIsFilterOpen(!isFilterOpen)}
          className="mb-4 px-4 py-2 bg-gray-200 text-gray-700 rounded-md hover:bg-gray-300"
        >
          {isFilterOpen ? "Hide Filters" : "Show Filters"}
        </button>

        {/* Collapsible Filter Menu (only Date Range) */}
        {isFilterOpen && (
          <div className="mb-6 space-y-2 p-4 border rounded-md bg-white">
            <div className="max-w-full sm:max-w-sm">
              <label className="font-medium text-gray-700 mb-1 block">Date Range:</label>
              <DatePicker
                selectsRange
                startDate={startDate}
                endDate={endDate}
                onChange={(update: [Date | null, Date | null]) => {
                  setStartDate(update[0]);
                  setEndDate(update[1]);
                }}
                isClearable={true}
                className="p-2 border border-gray-300 rounded-md w-full"
              />
            </div>
          </div>
        )}

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
                    <div className="flex justify-between items-start">
                      <div className="flex-1">
                        <p className="font-bold">{article.title}</p>
                        {previewArticleId === article.id && (
                          <div className="mt-1 text-sm text-gray-700 space-y-1">
                            <p><span className="font-semibold">Authors:</span> {article.authors}</p>
                            <p><span className="font-semibold">Journal:</span> {article.journal}</p>
                            <p><span className="font-semibold">Published:</span> {article.pub_date}</p>
                            <p><span className="font-semibold">Abstract:</span> {article.abstract}</p>
                            <a 
                              href={article.link} 
                              target="_blank" 
                              rel="noopener noreferrer"
                              className="text-blue-500 hover:underline"
                            >
                              View on PubMed
                            </a>
                          </div>
                        )}
                      </div>
                      <div className="flex flex-col items-end space-y-1">
                        <button
                          onClick={() =>
                            setPreviewArticleId(previewArticleId === article.id ? null : article.id)
                          }
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