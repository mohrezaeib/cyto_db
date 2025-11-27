import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react-swc';
import path from 'path';

export default defineConfig(({ mode }) => {
  // 👇 path to the repo root where .env lives
  const rootEnvDir = path.resolve(__dirname, '..');
  const env = loadEnv(mode, rootEnvDir, '');

  return {
    plugins: [react()],
    envPrefix: ['VITE_', 'CYTOCHAL_'],

    // 👇 tell Vite: "env files are one folder up"
    envDir: rootEnvDir,

    server: {
      host: '0.0.0.0',
      port: parseInt(env.CYTOCHAL_FRONTEND_PORT) || 3000,
    },
    build: {
      outDir: 'dist',
      sourcemap: false,
      assetsDir: 'assets',
    },
    base: env.CYTOCHAL_FRONTEND_BASE_URL || '/',
  };
});
