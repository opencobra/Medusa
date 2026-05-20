# Releasing medusa-cobra

This project publishes to [PyPI as `medusa-cobra`](https://pypi.org/project/medusa-cobra/).
Releases are cut from the `development` branch.

## Per-release steps

1. **Pick a version.** Compare against the current PyPI release and
   apply [SemVer](https://semver.org). New public API → minor bump.
   Bug fix only → patch bump.

2. **Update `CHANGELOG.md`.** Add a new section with the version and
   today's date. Group entries under `Added` / `Changed` / `Deprecated` /
   `Removed` / `Fixed` / `Infrastructure`.

3. **Update `setup.py`.** Bump the `version=` literal and the `vX.Y.Z`
   in `download_url`.

4. **Open a release PR** to `development` titled
   `Release X.Y.Z`. Get review, merge.

5. **Tag the merge commit and push.** On `development` after the merge:

   ```bash
   git pull --ff-only
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

   Pushing a `v*` tag triggers `.github/workflows/publish.yml`, which
   builds an sdist + wheel and uploads to PyPI via trusted publishing.

6. **Create a GitHub Release** from the tag. Paste the new
   `CHANGELOG.md` section as the body so users browsing the release page
   see what changed.

7. **Verify** the new version is live:

   ```bash
   pip install --upgrade medusa-cobra
   python -c "import medusa; print(medusa.__version__)"
   ```

## One-time setup: PyPI trusted publishing

The `publish.yml` workflow uses [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/)
(OIDC, no long-lived tokens). Before the first tag-triggered release
you must register the publisher on PyPI:

1. Log in to <https://pypi.org/manage/project/medusa-cobra/settings/publishing/>
   as a maintainer.
2. Add a **pending publisher** with:
   - Owner: `opencobra`
   - Repository: `medusa`
   - Workflow filename: `publish.yml`
   - Environment name: `pypi`
3. Save. The next tag push to `v*` will authenticate via OIDC.

## Manual fallback

If the GHA publish job is broken or unavailable, release manually from
a clean `development` checkout:

```bash
python -m pip install --upgrade build twine
python -m build
twine check dist/*
twine upload dist/*
git tag vX.Y.Z
git push origin vX.Y.Z
```
