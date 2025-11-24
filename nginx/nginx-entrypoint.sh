#!/bin/sh
set -e

# Defaults if env vars are not set
: "${BACKEND_PORT}"
: "${FRONTEND_PORT}"

echo "Starting Nginx Configuration..."
echo "  Frontend: http://frontend:${FRONTEND_PORT}"
echo "  Backend: http://backend:${BACKEND_PORT}"
echo "  Nginx Port: 80"

# Replace placeholders in template with env vars
envsubst '$BACKEND_PORT $FRONTEND_PORT' \
  < /etc/nginx/nginx.conf.template \
  > /etc/nginx/nginx.conf

echo "Nginx configuration completed successfully!"
