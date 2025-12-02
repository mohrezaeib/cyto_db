import os
from whoosh import index
from whoosh.fields import Schema, TEXT, ID, STORED
from whoosh.qparser import QueryParser
from whoosh.analysis import StemmingAnalyzer

class SearchEngine:
    def __init__(self, index_dir="search_index"):
        self.index_dir = index_dir
        self.schema = Schema(
            id=ID(stored=True, unique=True),
            Compound=TEXT(analyzer=StemmingAnalyzer(), stored=True),
            Origin=TEXT(analyzer=StemmingAnalyzer(), stored=True),
            Reference=TEXT(analyzer=StemmingAnalyzer(), stored=True),
            Smiles=TEXT(analyzer=StemmingAnalyzer(), stored=True),
            Subclass=TEXT(analyzer=StemmingAnalyzer(), stored=True),
            Molecular_Formula=TEXT(analyzer=StemmingAnalyzer(), stored=True)
        )
        self._setup_index()
    
    def _setup_index(self):
        """Create or open the search index"""
        if not os.path.exists(self.index_dir):
            os.makedirs(self.index_dir)
        
        if not index.exists_in(self.index_dir):
            self.ix = index.create_in(self.index_dir, self.schema)
        else:
            self.ix = index.open_dir(self.index_dir)
    
    def search(self, query_text, field_names=None):
        """Search across specified fields"""
        if not query_text:
            return []
        
        # Default fields to search
        if not field_names:
            field_names = ['Compound', 'Origin', 'Reference', 'Smiles', 'Subclass', 'Molecular_Formula']
        
        results = []
        with self.ix.searcher() as searcher:
            # Search each field and combine results
            for field in field_names:
                parser = QueryParser(field, schema=self.schema)
                query = parser.parse(query_text)
                field_results = searcher.search(query, limit=100)
                results.extend([hit['id'] for hit in field_results])
        
        unique_results = list(set(results))
        print(f"Search for '{query_text}': {len(unique_results)} results")
        return unique_results
    
    def index_document(self, doc_id, fields_dict):
        """Add or update a document in the index"""
        writer = self.ix.writer()
        
        doc_data = {
            'id': str(doc_id),
            'Compound': fields_dict.get('Compound', ''),
            'Origin': fields_dict.get('Origin', ''),
            'Reference': fields_dict.get('Reference', ''),
            'Smiles': fields_dict.get('Smiles', ''),
            'Subclass': fields_dict.get('Subclass', ''),
            'Molecular_Formula': fields_dict.get('Molecular Formula', '')
        }
        writer.update_document(**doc_data)
        writer.commit()

# Global instance
search_engine = SearchEngine()