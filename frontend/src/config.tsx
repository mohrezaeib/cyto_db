// config.tsx

interface Config {
  apiUrl: string;
  isDevelopment: boolean;   // Vite flag
  apiHost: string;
  apiPort: string;
  apiBasePath: string;
  staticBaseUrl: string;
  microscopyImageBaseUrl: string;
  moleculeImageBaseUrl: string;
}

// helper to read env vars safely
const getEnvVar = (key: string, defaultValue = ''): string => {
  const value = (import.meta.env as any)[key];
  return typeof value === 'string' && value.length > 0 ? value : defaultValue;
};

const getConfig = (): Config => {
  const isDevelopment = import.meta.env.DEV; // Vite built-in

  const apiHost = getEnvVar('CYTOCHAL_API_HOST', 'localhost');
  const apiPort = getEnvVar('CYTOCHAL_BACKEND_PORT', '5000');
  const apiBasePathRaw = getEnvVar('CYTOCHAL_API_BASE_PATH', '/cytochal-api');

  // normalize base path
  const apiBasePath = apiBasePathRaw.startsWith('/')
    ? apiBasePathRaw
    : `/${apiBasePathRaw}`;

  // from .env: CYTOCHAL_IS_DEV=true/false
  const isDevEnv =
    getEnvVar('CYTOCHAL_IS_DEV', 'false').toLowerCase() === 'true';

  console.log('CYTOCHAL_IS_DEV =', getEnvVar('CYTOCHAL_IS_DEV', 'false'));

  // API base (you already saw this working)
  const apiUrl = isDevEnv
    ? `http://${apiHost}:${apiPort}${apiBasePath}` // dev → backend:5000
    : apiBasePath;                                // prod → nginx relative

  // STATIC files base
  const staticBaseUrl = isDevEnv
    ? `http://${apiHost}:${apiPort}/static`       // dev → backend:5000/static
    : `/static`;                                  // prod → nginx /static

  // Specific image bases
  const microscopyImageBaseUrl = `${staticBaseUrl}/microscopy_images/`;
  const moleculeImageBaseUrl = `${staticBaseUrl}/molecule_structures_image/`;

  return {
    apiUrl,
    isDevelopment,
    apiHost,
    apiPort,
    apiBasePath,
    staticBaseUrl,
    microscopyImageBaseUrl,
    moleculeImageBaseUrl,
  };
};

export const config = getConfig();
