#!/bin/bash

# Define variables
BRANCH_NAME="fair-submission"
COMMIT_MSG="Initial Release: Alzheimer's Tau Protein Docking Simulation v1.0"

echo "Starting Git History Reset to Orphan Branch: $BRANCH_NAME"

# Check if we are in a git repo
if [ ! -d ".git" ]; then
    echo "Error: Not a git repository."
    exit 1
fi

# Create orphan branch
git checkout --orphan "$BRANCH_NAME"

if [ $? -ne 0 ]; then
    echo "Error creating orphan branch."
    exit 1
fi

echo "Orphan branch created."

# Add all files
git add .

echo "Files added."

# Commit
git commit -m "$COMMIT_MSG"

if [ $? -ne 0 ]; then
    echo "Error committing files."
    exit 1
fi

echo "Commit successful."

echo "--------------------------------------------------------"
echo "Git history reset complete on branch '$BRANCH_NAME'."
echo "To push this new clean history to your remote repository, run:"
echo ""
echo "  git push -f origin $BRANCH_NAME"
echo ""
echo "OR to overwrite main (DANGEROUS - deletes old history on remote main):"
echo ""
echo "  git checkout $BRANCH_NAME"
echo "  git branch -D main"
echo "  git branch -m main"
echo "  git push -f origin main"
echo "--------------------------------------------------------"
