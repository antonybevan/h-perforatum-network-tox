import subprocess
import sys

def run_step(command, step_name):
    print(f"--- Running {step_name} ---")
    try:
        subprocess.check_call(command, shell=True)
        print(f"--- {step_name} Passed ---\n")
    except subprocess.CalledProcessError:
        print(f"--- {step_name} FAILED ---\n")
        sys.exit(1)

def main():
    print("Starting CI Build...")
    
    # 1. Install/Check Dependencies (Optional in CI but good for local/agent consistency)
    # run_step(f"{sys.executable} -m pip install -r requirements.txt", "Dependency Check")

    # 2. Linting (Industry Standard Style Enforcement)
    # Checks for bugs, syntax errors, and style violations
    run_step(f"{sys.executable} -m ruff check .", "Linting")

    # 3. Run Tests
    # Using python -m pytest to ensure it uses the current python environment
    run_step(f"{sys.executable} -m pytest", "Unit Tests")

    print("Build Successful!")

if __name__ == "__main__":
    main()
