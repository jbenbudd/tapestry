// frontend/src/App.tsx
import React, { useState } from 'react';
import SearchSection from './components/SearchSection';
import BuildSection from './components/BuildSection';
import ChatSection from './components/ChatSection';

export interface Article {
  id: string;
  title: string;
  abstract: string;
  authors?: string;
  journal?: string;
  pub_date?: string;
  link?: string;
}

const App: React.FC = () => {
  const [selectedArticles, setSelectedArticles] = useState<Article[]>([]);

  return (
    <div className="min-h-screen bg-gray-100 p-4">
      <h1 className="text-3xl font-bold text-center mb-6">
        Tapestry: A retrieval-augmented generation engine for academic research
      </h1>
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4 md:h-[80vh]">
        <SearchSection 
          selectedArticles={selectedArticles} 
          setSelectedArticles={setSelectedArticles}
        />
        <BuildSection 
          selectedArticles={selectedArticles} 
          setSelectedArticles={setSelectedArticles} 
        />
        <ChatSection />
      </div>
    </div>
  );
};

export default App;