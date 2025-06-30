"""
GWAS QC Pipeline Utilities
Author: Ren F (tt7)
Copyright (c) 2025
"""

#!/usr/bin/env python3

import os
import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict, Any
import json
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_dependency(tool: str) -> bool:
    """Check if a required tool is installed."""
    try:
        subprocess.run(['which', tool], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        logger.error(f"Required tool '{tool}' is not installed")
        return False

def validate_file_exists(file_path: str, file_type: str) -> bool:
    """Validate that a file exists and is readable."""
    path = Path(file_path)
    if not path.exists():
        logger.error(f"{file_type} file not found: {file_path}")
        return False
    if not os.access(path, os.R_OK):
        logger.error(f"{file_type} file is not readable: {file_path}")
        return False
    return True

def validate_vcf_format(vcf_file: str) -> bool:
    """Validate VCF file format using bcftools."""
    try:
        result = subprocess.run(
            ['bcftools', 'view', '-h', vcf_file],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            logger.error(f"Invalid VCF file format: {result.stderr}")
            return False
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error validating VCF file: {e}")
        return False

def safe_json_dump(data: Dict[str, Any], output_file: str) -> bool:
    """Safely write data to a JSON file."""
    try:
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=4)
        return True
    except Exception as e:
        logger.error(f"Error writing to JSON file: {e}")
        return False

def get_file_size(file_path: str) -> Optional[int]:
    """Get file size in bytes."""
    try:
        return os.path.getsize(file_path)
    except OSError as e:
        logger.error(f"Error getting file size: {e}")
        return None

def check_disk_space(path: str, required_space_mb: int) -> bool:
    """Check if there's enough disk space available."""
    try:
        stat = os.statvfs(path)
        free_space = stat.f_frsize * stat.f_bavail
        required_space = required_space_mb * 1024 * 1024
        return free_space >= required_space
    except OSError as e:
        logger.error(f"Error checking disk space: {e}")
        return False

def safe_save_df(df, file_path, index=False, sep='\t'):
    """Safely save a DataFrame to a file, creating the directory if needed."""
    try:
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        df.to_csv(file_path, index=index, sep=sep)
        logger.info(f"Successfully saved DataFrame to {file_path}")
    except Exception as e:
        logger.error(f"Failed to save DataFrame to {file_path}: {e}")
        raise

def safe_load_df(file_path, sep='\t'):
    """Safely load a DataFrame from a file."""
    try:
        validate_file_exists(file_path, "input DataFrame")
        df = pd.read_csv(file_path, sep=sep)
        logger.info(f"Successfully loaded DataFrame from {file_path}")
        return df
    except Exception as e:
        logger.error(f"Failed to load DataFrame from {file_path}: {e}")
        raise 