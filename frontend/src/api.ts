// api.ts
import axios from 'axios';
import { FilterParams, Compound, ApiResponse } from './types';
import { config } from './config';

// Use config that works in both dev and prod
const API_BASE_URL = config.apiUrl;
console.log({API_BASE_URL});

export const fetchCompounds = async (params: FilterParams = {}): Promise<ApiResponse> => {
  try {
    const response = await axios.get<ApiResponse>(`${API_BASE_URL}/items`, {
      params: {
        page: 1,
        per_page: 50,
        ...params
      }
    });
    return response.data;
  } catch (error) {
    console.error('Error fetching compounds:', error);
    throw error;
  }
};