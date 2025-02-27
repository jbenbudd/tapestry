import React, { useState } from 'react';
import { Article } from '../App';

interface BuildSectionProps {
  selectedArticles: Article[];
  setSelectedArticles: React.Dispatch<React.SetStateAction<Article[]>>;
}

const BuildSection: React.FC<BuildSectionProps> = ({ selectedArticles, setSelectedArticles }) => {
  const [isBuilding, setIsBuilding] = useState<boolean>(false);
  const [progress, setProgress] = useState<number>(0);

  const handleBuildDatabase = () => {
    setIsBuilding(true);
    setProgress(0);
    // Simulate progress with an interval
    const interval = setInterval(() => {
      setProgress(prev => {
        if (prev >= 100) {
          clearInterval(interval);
          setIsBuilding(false);
          return 100;
        }
        return prev + 10;
      });
    }, 500);
  };

  const removeArticle = (id: string) => {
    setSelectedArticles(selectedArticles.filter(article => article.id !== id));
  };

  return (
    <div className="bg-white p-4 rounded-md shadow-md flex flex-col h-[50vh] md:h-full min-h-0">
      <div>
        <h2 className="text-xl font-semibold mb-4">Build</h2>
        <p className="mb-2">
          You have {selectedArticles.length} document(s) selected.
        </p>
        <button
          onClick={handleBuildDatabase}
          className="mb-4 px-4 py-2 bg-green-500 text-white rounded-md hover:bg-green-600"
        >
          Build Database
        </button>
        {isBuilding && (
          <div className="w-full bg-gray-300 rounded-full h-4">
            <div
              className="bg-green-500 h-4 rounded-full"
              style={{ width: `${progress}%` }}
            />
          </div>
        )}
      </div>
      <div className="mt-4 overflow-y-auto flex-1 min-h-0 pr-4">
        <h3 className="text-lg font-semibold mb-2">Selected Articles</h3>
        {selectedArticles.length === 0 ? (
          <p>No articles added yet.</p>
        ) : (
          <ul>
            {selectedArticles.map(article => (
              <li key={article.id} className="border-b border-gray-300 py-2">
                <div className="flex justify-between items-center">
                  <p className="font-bold">{article.title}</p>
                  <button
                    onClick={() => removeArticle(article.id)}
                    className="text-sm text-red-500 hover:underline"
                  >
                    Remove
                  </button>
                </div>
                <div className="mt-1 text-sm text-gray-700">
                  <p>
                    <span className="font-semibold">Authors:</span> {article.authors}
                  </p>
                  <p>
                    <span className="font-semibold">Journal:</span> {article.journal}
                  </p>
                  <p>
                    <span className="font-semibold">Published:</span> {article.pub_date}
                  </p>
                  <a
                    href={article.link}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-500 hover:underline"
                  >
                    View on PubMed
                  </a>
                </div>
              </li>
            ))}
          </ul>
        )}
      </div>
    </div>
  );
};

export default BuildSection;