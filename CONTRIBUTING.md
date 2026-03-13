# Contributing to HERA

First off, thank you for considering contributing to HERA! It's people like you
that make the open source community such an amazing place to learn, inspire, and
create. Any contribution you make is **greatly appreciated**.

This document provides guidelines for contributing to the project.

## Collaborator Access

To maintain a clear overview of all ongoing developments and ensure efficient coordination,
**contributions are currently managed via collaborator access**.

Everyone is still warmly invited to contribute! If you would like to help improve HERA,
please reach out by opening a [Discussion](https://github.com/lerdmann1601/HERA-Matlab/discussions)
or an Issue to be added as a collaborator. This helps me keep track of who is working on what
and avoid duplicate efforts.

## Code of Conduct

This project and everyone participating in it is governed by the [Contributor
Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
By participating, you are expected to uphold this code. Please report
unacceptable behavior to the project maintainer.

## Getting Started

### Have a Question?

If you have questions about how to use HERA, please check the `README.md` and
existing documentation first. If you're still stuck, feel free to open a
[Discussion](https://github.com/lerdmann1601/HERA-Matlab/discussions) or an
Issue with the label `question`.

### Found a Bug?

If you spot a bug, please help us by reporting it!

1. **Search existing issues** to see if it has already been reported.
2. If not, open a new Issue using the **Bug Report template**.
3. Provide all requested details in the template:
    * Your Operating System and MATLAB version.
    * Steps to reproduce the error.
    * Screenshots or code snippets if applicable.

## How to Contribute

### 1. Fork and Clone

Fork the repository to your own GitHub account and then clone it to your local machine:

```bash
git clone https://github.com/YOUR-USERNAME/HERA-Matlab.git
```

### 2. Create a Branch

Create a new branch for your feature or fix. This keeps your changes organized and separate from the main codebase.

```bash
git checkout -b feature/AmazingFeature
```

### 3. Make Your Changes

Implement your feature or bug fix.

* **Code Style**: Please follow standard MATLAB and Python coding conventions.
* **Documentation**: Update the `README.md` or add comments if you are introducing new functionality.

### 4. Run Tests

Before submitting, please ensure that your changes don't break existing functionality. Run the included test suite:

```matlab
run('tests/run_HERA_tests.m')
```

### 5. Commit and Push

Commit your changes with a clear and descriptive message:

```bash
git commit -m "Add some AmazingFeature"
git push origin feature/AmazingFeature
```

### 6. Open a Pull Request

Go to the original HERA repository and open a Pull Request.

* **Use the Pull Request template**: GitHub will automatically populate the description field with our template.
Please fill it out completely.
* Describe your changes clearly.
* Link to any relevant issues (e.g., "Fixes #123").
* Wait for review! I will do my best to review your contribution as soon as possible.

---

**Thank you for your hard work and support!**
