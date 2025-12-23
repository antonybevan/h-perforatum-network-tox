# Git Repository Setup Instructions

## ✅ Local Repository Ready

**Initial commit created:** `cfe5c8f`  
**Files committed:** 47  
**Repository:** e:\network_pharmacology\h-perforatum-net-tox

---

## Next Steps: Push to GitHub

### Option 1: Using GitHub Website (Recommended)

1. **Create Repository on GitHub:**
   - Go to https://github.com/new
   - Repository name: `h-perforatum-network-tox`
   - Description: "Network pharmacology analysis demonstrating Hyperforin's 26x higher per-target hepatotoxic influence"
   - Choose: Public or Private
   - **DO NOT** initialize with README, .gitignore, or license (we already have these)
   - Click "Create repository"

2. **Push to GitHub:**
   ```bash
   cd e:\network_pharmacology\h-perforatum-net-tox
   
   # Add remote (replace antonybevan)
   git remote add origin https://github.com/antonybevan/h-perforatum-network-tox.git
   # Push to GitHub
   git branch -M main
   git push -u origin main
   ```

### Option 2: Using GitHub CLI

```bash
# Install GitHub CLI if not installed
# Then:
gh repo create h-perforatum-network-tox --public --source=. --remote=origin
git push -u origin main
```

---

## Recommended Repository Settings

### Topics (Add these on GitHub)
```
network-pharmacology
drug-toxicity
systems-biology
bioinformatics
python
networkx
```

### Repository Description
```
Network pharmacology analysis of H. perforatum hepatotoxicity using RWR and tissue-specific networks. Demonstrates 26x higher per-target influence of Hyperforin vs Quercetin.
```

### .gitattributes (Optional - Add if needed)
```bash
# Handle line endings
* text=auto
*.py text
*.md text
*.csv text
*.txt text
```

---

## Post-Push Checklist

- [ ] Verify all files pushed correctly
- [ ] Check README displays properly
- [ ] Add repository topics
- [ ] Enable GitHub Issues (optional)
- [ ] Add CITATION.cff (optional, for academic citation)
- [ ] Set up branch protection rules (optional)

---

## Git Commands Reference

```bash
# Check status
git status

# View commit history
git log --oneline

# Create new branch
git checkout -b feature-name

# Push changes
git add .
git commit -m "Description"
git push
```

---

**Current Status:** Ready to push! 🚀
