/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly CYTOCHAL_API_HOST: string
  readonly CYTOCHAL_BACKEND_PORT: string
  readonly CYTOCHAL_API_BASE_PATH: string

  readonly CYTOCHAL_IS_DEV: string

  readonly VITE_API_BASE_URL: string
  readonly  staticBaseUrl: string;
 readonly microscopyImageBaseUrl: string;
 readonly moleculeImageBaseUrl: string;

}

interface ImportMeta {
  readonly env: ImportMetaEnv
}
