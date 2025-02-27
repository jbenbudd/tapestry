import React, { useState } from 'react';

const ChatSection: React.FC = () => {
  const [messages, setMessages] = useState<string[]>([]);
  const [input, setInput] = useState<string>('');
  const [selectedModel, setSelectedModel] = useState<string>('default-model');

  const handleSend = () => {
    if (!input.trim()) return;
    // Simulate sending a message and receiving an echo response
    setMessages([...messages, input, `Echo: ${input}`]);
    setInput('');
  };

  return (
    <div className="bg-white p-4 rounded-md shadow-md flex flex-col h-[50vh] md:h-full min-h-0">
      <h2 className="text-xl font-semibold mb-4">Chat</h2>
      <div className="mb-4">
        <label className="block text-gray-700 mb-1">Select Model:</label>
        <select
          value={selectedModel}
          onChange={(e) => setSelectedModel(e.target.value)}
          className="w-full p-2 border border-gray-300 rounded-md"
        >
          <option value="default-model">Default Model</option>
          <option value="model-a">Model A</option>
          <option value="model-b">Model B</option>
        </select>
      </div>
      <div className="flex-1 overflow-y-auto mb-4 border p-2 rounded-md min-h-0 pr-4">
        {messages.length === 0 ? (
          <p className="text-gray-500">No messages yet.</p>
        ) : (
          messages.map((msg, index) => (
            <p key={index} className="mb-2">{msg}</p>
          ))
        )}
      </div>
      <div className="flex">
        <input
          type="text"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          placeholder="Type your message..."
          className="flex-1 p-2 border border-gray-300 rounded-l-md"
        />
        <button
          onClick={handleSend}
          className="px-4 py-2 bg-blue-500 text-white rounded-r-md hover:bg-blue-600"
        >
          Send
        </button>
      </div>
    </div>
  );
};

export default ChatSection;