import React, { useState } from 'react';
import DatePicker from 'react-datepicker';
import 'react-datepicker/dist/react-datepicker.css';
import { Article } from '../App';

interface SearchSectionProps {
  selectedArticles: Article[];
  setSelectedArticles: React.Dispatch<React.SetStateAction<Article[]>>;
}

const SearchSection: React.FC<SearchSectionProps> = ({ selectedArticles, setSelectedArticles }) => {
  const [query, setQuery] = useState<string>('');
  const [startDate, setStartDate] = useState<Date | null>(null);
  const [endDate, setEndDate] = useState<Date | null>(null);
  const [results, setResults] = useState<Article[]>([]);
  const [previewArticleId, setPreviewArticleId] = useState<string | null>(null);
  const [isFilterOpen, setIsFilterOpen] = useState<boolean>(false);

  const handleSearch = async () => {
    if (!query) return;
    try {
      const formattedStartDate = startDate ? startDate.toISOString().split('T')[0] : '';
      const formattedEndDate = endDate ? endDate.toISOString().split('T')[0] : '';
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

  const handleClearResults = () => {
    setResults([]);
    setPreviewArticleId(null);
  };

  const handleAdd = (article: Article) => {
    if (!selectedArticles.find(a => a.id === article.id)) {
      setSelectedArticles([...selectedArticles, article]);
      setResults(results.filter(a => a.id !== article.id));
      if (previewArticleId === article.id) {
        setPreviewArticleId(null);
      }
    }
  };

  return (
    <div className="bg-white p-4 rounded-md shadow-md flex flex-col h-[50vh] md:h-full min-h-0">
      <h2 className="text-xl font-semibold mb-4">Search</h2>
      <div className="mb-4">
        <input
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Search PubMed articles..."
          className="w-full p-2 border border-gray-300 rounded-md"
        />
        <div className="flex space-x-2 mt-2">
          <button
            onClick={handleSearch}
            className="px-4 py-2 bg-blue-500 text-white rounded-md hover:bg-blue-600"
          >
            Search
          </button>
          <button
            onClick={handleClearResults}
            className="px-4 py-2 bg-red-500 text-white rounded-md hover:bg-red-600"
          >
            Clear
          </button>
        </div>
      </div>
      <button
        onClick={() => setIsFilterOpen(!isFilterOpen)}
        className="mb-4 px-4 py-2 bg-gray-200 text-gray-700 rounded-md hover:bg-gray-300"
      >
        {isFilterOpen ? 'Hide Filters' : 'Show Filters'}
      </button>
      {isFilterOpen && (
        <div className="mb-4">
          <label className="block text-gray-700 mb-1">Date Range:</label>
          <DatePicker
            selectsRange
            startDate={startDate}
            endDate={endDate}
            onChange={(update: [Date | null, Date | null]) => {
              setStartDate(update[0]);
              setEndDate(update[1]);
            }}
            isClearable={true}
            className="w-full p-2 border border-gray-300 rounded-md"
          />
        </div>
      )}
      <div className="mt-4 overflow-y-auto flex-1 min-h-0 pr-4">
        <h3 className="text-lg font-semibold">Results</h3>
        {results.length === 0 ? (
          <p>No results yet.</p>
        ) : (
          <ul>
            {results.map(article => (
              <li key={article.id} className="border-b border-gray-200 py-2">
                <div className="flex justify-between items-center">
                  <p className="font-bold">{article.title}</p>
                  <div>
                    <button
                      onClick={() =>
                        setPreviewArticleId(
                          previewArticleId === article.id ? null : article.id
                        )
                      }
                      className="text-sm text-blue-500 hover:underline mr-2"
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
                {previewArticleId === article.id && (
                  <div className="mt-1 text-sm text-gray-700">
                    <p><span className="font-semibold">Authors:</span> {article.authors}</p>
                    <p><span className="font-semibold">Journal:</span> {article.journal}</p>
                    <p><span className="font-semibold">Published:</span> {article.pub_date}</p>
                    <p>{article.abstract}</p>
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
              </li>
            ))}
          </ul>
        )}
      </div>
    </div>
  );
};

export default SearchSection;