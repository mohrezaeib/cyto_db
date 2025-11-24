#!/bin/bash

# Run the SDF parser first to generate required JSON data
echo "Running SDF parser to generate database files..."
python cyto_db_shiny_app/parse_sdf.py

# Check if parser ran successfully
if [ $? -eq 0 ]; then
    echo "SDF parser completed successfully"
    
    # Wait a moment to ensure file is fully written
    sleep 2
    
    # Verify the data file was created
    if [ -f "static/data.json" ]; then
        echo "Data file created successfully at static/data.json"
        echo "File size: $(stat -c%s static/data.json) bytes"
    else
        echo "ERROR: Data file not found!" >&2
        exit 1
    fi
else
    echo "Error: SDF parser failed" >&2
    exit 1
fi

# Start the Flask application
WORKERS=$(python -c "import multiprocessing as m; print(2 * m.cpu_count() + 1)")
echo "Starting Flask application with Gunicorn..."
echo "Using $WORKERS Gunicorn workers..."
exec gunicorn --bind 0.0.0.0:5000 --workers "$WORKERS" --preload app:app