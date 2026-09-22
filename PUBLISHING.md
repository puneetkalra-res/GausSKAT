# Publishing this directory on GitHub

Create an empty GitHub repository named `GausSKAT`. From this directory, run:

```bash
git init -b main
git add .
git commit -m "Initial independent GausSKAT implementation"
git remote add origin https://github.com/puneetkalra-res/GausSKAT.git
git push -u origin main
```

The included GitHub Actions workflow will run `R CMD check` after the first
push. Review that check before sharing the repository with the editor.

No new GitHub account is technically required. Use a separate account or a
private repository only if journal blinding or your preferred separation from
personal repositories requires it.
