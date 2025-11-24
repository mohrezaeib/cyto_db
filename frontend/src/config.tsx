// config.ts
interface Config {
  apiUrl: string;
  isDevelopment: boolean;
  apiHost: string;
  apiPort: string;
  apiBasePath: string;
}

// Type-safe environment variable access
const getEnvVar = (key: string, defaultValue: string = ''): string => {
  const value = import.meta.env[key];
  return value ? String(value) : defaultValue;
};

const getConfig = (): Config => {
  const isDevelopment = import.meta.env.DEV;
  
  const apiHost = getEnvVar('VITE_CYTOCHAL_API_HOST', 'localhost');
  const apiPort = getEnvVar('VITE_CYTOCHAL_BACKEND_PORT', '5000');
  const apiBasePath = getEnvVar('VITE_CYTOCHAL_API_BASE_PATH', '/cytochal-api');

  if (isDevelopment) {
    // Development mode - use proxy
    const apiUrl = `${apiBasePath}`; // Proxy will handle the host/port
    return {
      apiUrl,
      isDevelopment: true,
      apiHost,
      apiPort,
      apiBasePath
    };
  } else {
    // Production mode - use full URL from environment variables
    const apiUrl = getEnvVar('VITE_CYTOCHAL_API_FULL_URL') || 
                   `http://${apiHost}:${apiPort}${apiBasePath}`;
    return {
      apiUrl,
      isDevelopment: false,
      apiHost,
      apiPort,
      apiBasePath
    };
  }
};

export const config = getConfig();

// Helper function to log config (useful for debugging)
export const logConfig = () => {
  if (config.isDevelopment) {
    console.log('App Configuration:', {
      mode: 'development',
      apiUrl: config.apiUrl,
      apiHost: config.apiHost,
      apiPort: config.apiPort,
      apiBasePath: config.apiBasePath
    });
  }
};