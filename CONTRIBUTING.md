# CONTRIBUTING to variantCell

Thank you for your interest in contributing to variantCell! This document outlines how you can contribute to this R package for single-cell SNP analysis.

## Table of Contents

1.  [Code of Conduct](#code-of-conduct)
2.  [Getting Started](#getting-started)
3.  [Development Workflow](#development-workflow)
4.  [Pull Request Process](#pull-request-process)
5.  [Style Guidelines](#style-guidelines)
6.  [Testing](#testing)
7.  [Issue Reporting](#issue-reporting)
8.  [Feature Requests](#feature-requests)
9.  [Contact](#contact)

## Code of Conduct {#code-of-conduct}

By participating in this project, you agree to uphold our Code of Conduct (we follow the [Contributor Covenant](https://www.contributor-covenant.org/)). Please ensure all interactions are respectful, inclusive, and professional.

## Getting Started {#getting-started}

### Prerequisites

-   R (latest stable version recommended)
-   Core dependencies: R6, data.table, Matrix, ggplot2, cowplot, GenomicRanges, IRanges, AnnotationHub, matrixStats, ensembldb, circlize, ComplexHeatmap, parallel, doParallel, foreach
-   Optional: Seurat (for Seurat objects), SingleCellExperiment (for SCE objects)

### Fork and Clone

1.  Fork the repository on GitHub
2.  Clone your fork locally:

``` bash
git clone https://github.com/YOUR-USERNAME/variantCell.git
cd variantCell
```

3.  Add the original repository as an upstream remote:

```         
bash
git remote add upstream https://github.com/potterae/variantCell.git
```

4.  Install the package in development mode:

```         
bash
rdevtools::install_local("path/to/variantCell", force = TRUE)
```

### Development Workflow {#development-workflow}

1.  Create a new branch for your feature or bugfix:

```         
bash
git checkout -b feature/your-feature-name
```

or

```         
bash
git checkout -b fix/issue-you-are-fixing
```

2.  Make your changes, following the style guidelines below Test your changes with example data (see Testing)

3.  Test your changes

4.  Commit your changes with descriptive commit messages

```         
bash
git commit -m "Add functionality to analyze donor-specific SNPs"
```

5.  Push to your fork:

```         
bash
git push origin feature/your-feature-name
```

6.  Submit a pull request (see Pull Request Process)

### Pull Request Process {#pull-request-process}

1.  Update the README.md and documentation with details of changes if appropriate
2.  Update the version number in relevant files following semantic versioning
3.  Submit your pull request to the main repository
4.  Respond to any feedback or questions during the review process
5.  After approval, a maintainer will merge your PR

### Style Guidelines {#style-guidelines}

We follow the tidyverse style guide for R code with the following specifics:

1.  Use camelCase for function names within R6 classes, following the existing pattern
2.  Add comprehensive Roxygen documentation for all public functions
3.  Include examples in the documentation
4.  For R6 class methods, follow the documentation format shown in the existing code
5.  Indent with 2 spaces
6.  Maximum line length of 100 characters
7.  Include relevant error handling and input validation

### Testing {#testing}

1.  Create test cases for new functionality
2.  Verify that your code works with example data
3.  If adding new features, provide a small example dataset if possible
4.  Test with both transplant and non-transplant modes if applicable
5.  For performance-critical functions, include benchmarking

Issue Reporting {#issue-reporting}

When reporting issues, please include:

1.  A clear description of the issue
2.  Steps to reproduce the problem
3.  Expected behavior
4.  Actual behavior
5.  R version, package version, and OS information
6.  If possible, a minimal reproducible example

Use the issue templates provided in the repository.

###Feature Requests {#feature-requests}

Feature requests are welcome! Please provide:

1.  A clear description of the proposed feature
2.  The motivation or use case for this feature
3.  If possible, a sketch or example of how it might work
4.  Whether you're interested in implementing it yourself

Contact {#contact}

For questions about contributing, please:

Open an issue on GitHub Contact the package maintainer: Andrew Potter

Thank you for your interest in improving variantCell!
