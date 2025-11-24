import { defineConfig, loadEnv } from 'vite'
import react from '@vitejs/plugin-react-swc'

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '')
  
  return {
    plugins: [react()],
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
  }
})