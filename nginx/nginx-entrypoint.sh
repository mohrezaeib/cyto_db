#!/bin/sh
set -e

# Defaults if env vars are not set
: "${BACKEND_PORT}"
: "${FRONTEND_PORT}"
: "${SERVER_NAMES:=localhost}" # Default to localhost if not set


echo "Starting Nginx Configuration..."
echo "  Frontend: http://frontend:${FRONTEND_PORT}"
echo "  Backend: http://backend:${BACKEND_PORT}"
echo "  Server Names: ${SERVER_NAMES}"
echo "  Nginx Port: 80"

# Replace placeholders in template with env vars
envsubst '$BACKEND_PORT $FRONTEND_PORT' \
  < /etc/nginx/nginx.conf.template \
  > /etc/nginx/nginx.conf

echo "Nginx configuration completed successfully!"
