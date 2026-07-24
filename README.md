## configuration of toolshed galaxy packages

Galaxy tool wrappers published to the Galaxy Tool Shed. Each top-level
directory is one Tool Shed repository (see `CLAUDE.md` for details).

### Publishing to the main Tool Shed

Done manually by the maintainer. From inside the tool's directory, with the
Tool Shed API key in `$KEY`:

```
planemo shed_update --shed_target toolshed --shed_key $KEY --owner petr-novak .
```

Test on the testtoolshed sandbox first (`--shed_target testtoolshed
--owner petrn`); only push to the main toolshed once `planemo test` passes.
