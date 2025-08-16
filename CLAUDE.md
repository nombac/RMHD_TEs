# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains analysis tools for ZEUS radiation hydrodynamics simulations, supporting both IDL and Python workflows. The codebase focuses on analyzing accretion disk physics, thermal instabilities, and magnetorotational turbulence through S-curve analysis and other diagnostic tools.

## Development Commands

### Fortran Simulation (test/src/)
- Build simulation: `cd test/src && make SYSTYPE=macosx TARG=<target_name>`
- Clean build files: `cd test/src && make clean`
- Available system types: `macosx`, `xt`, `mstatx00`, `altix`, `teragrid`

### Python Analysis
- Run S-curve analysis: `python s_curve_all.py`
- Required packages: `pip install numpy matplotlib`

### IDL Analysis
- Initialize IDL environment: Load `startup.pro` (defines physical constants and plotting parameters)
- Main analysis procedures: `have.pro`, `prof.pro`, `vave.pro`, `slice.pro`

## Code Architecture

### Core Components

**Simulation Code (test/src/)**
- Fortran 90 modules for radiation magnetohydrodynamics
- MPI-parallel with radiation transport and magnetic fields
- Modular design: `field_module`, `radiation_module`, `transport_module`, etc.
- Preprocessor macros control physics: COMPTON, MPI, SHEAR, DAMPING

**IDL Analysis Tools**
- `zeus_param.pro`: Extracts simulation parameters from Makefiles and input files
- `startup.pro`: Defines physical constants and plotting symbols as system variables
- `have.pro`: Horizontally-averaged diagnostic analysis
- `prof.pro`: Radial profile analysis
- `slice.pro`: 2D slice visualization
- Data readers: `read_dump.pro`, `readu.pro`, `read_qes.pro`

**Python Analysis (newer)**
- `s_curve_all.py`: Main S-curve plotting (converted from IDL)
- `zeus_utils.py`: Binary data reading utilities
- `target_config.py`: Simulation target definitions
- `constants.py`: Physical constants matching IDL definitions

### Data Flow
1. Fortran simulation generates binary dumps in `data/` directory
2. IDL/Python tools read binary data using custom readers
3. `zeus_param.pro` extracts metadata from Makefile and input files
4. Analysis scripts process time series and spatial data
5. Output plots saved to `eps/` or specified directories

### Key Conventions

**File Organization**
- Simulation source: `test/src/`
- Analysis scripts: Root directory
- Historical data: `old/` subdirectory
- Specialized analysis: `mnras.2017/` subdirectory

**IDL Conventions**
- Use `COMPILE_OPT IDL2` for modern IDL syntax
- Physical constants defined as system variables (e.g., `!MSOL`, `!SIGMAB`)
- Data paths constructed relative to simulation target directories
- Binary data typically big-endian format

**Python Conventions**
- Modern Python 3 with numpy/matplotlib
- Binary data reading uses struct module with big-endian format
- Configuration separated into dedicated modules
- Output as PDF instead of PostScript

### Data Formats
- Simulation dumps: Custom binary format (Fortran unformatted)
- Parameter files: ASCII text (`z3dinput`, `ipara.data`, etc.)
- Time series: Binary arrays (`resolh.data`, `vave.data`)
- Endianness: Big-endian for cross-platform compatibility

## Physics Context

This codebase analyzes thermal instability in accretion disks through:
- S-curve analysis (surface density vs. effective temperature)
- Alpha parameter evolution (turbulent stress analysis)  
- Vertical averaging and radial profiles
- Time series analysis of thermodynamic quantities
- Radiation cooling curve diagnostics

The `test/` directory contains a complete radiation MHD simulation setup with job scripts for various HPC systems.

## Language Preference

The user prefers to communicate in Japanese. Use Japanese for all interactions and responses when working with this codebase.

## Work History and Logs

Check the `logs/` directory for detailed work logs that document previous development sessions. These logs contain:
- IDL to Python conversion progress
- File dependencies and implementation details  
- Test results and validation steps
- Remaining tasks and next steps

When starting work, review the most recent log file to understand the current state and continue from where previous work left off. This helps maintain continuity across different Claude Code sessions.

## Git Commit Guidelines

When committing changes to GitHub:
- Do NOT include Claude as co-author in commit messages
- Remove the "Co-Authored-By: Claude <noreply@anthropic.com>" line from commit messages
- Use standard commit message format without AI attribution to avoid potential copyright complications
- CLAUDE.md should NOT be included in the GitHub repository

## Logging Guidelines

When the user requests "ログを記録して" (record logs):
- Document all actions performed in the current Claude Code session
- Write the log in a way that is easy to understand when reviewed later
- Follow the granularity, format, and file naming conventions of existing logs in the `logs/` directory
- Use file naming convention: `YYYY-MMDD-hhmm.md` (e.g., `2025-0815-1040.md`)
- If the user explains what they want to do next after requesting logs, add a "## 次にやりたいこと" section at the end of the log with their explanation. If nothing is mentioned, leave this section empty.
- Report completion once the log has been written