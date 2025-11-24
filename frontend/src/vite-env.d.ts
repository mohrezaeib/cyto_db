/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_CYTOCHAL_API_HOST: string
  readonly VITE_CYTOCHAL_BACKEND_PORT: string
  readonly VITE_CYTOCHAL_API_BASE_PATH: string
  readonly VITE_CYTOCHAL_API_FULL_URL: string

}

interface ImportMeta {
  readonly env: ImportMetaEnv
}