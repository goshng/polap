# 2117 ctags -R --languages=Bash --kinds-Bash=f -f .tags .
# 2118 ctags --version | head -1
# 2119 ctags -R --languages=sh --kinds-sh=f -f .tags .
# 2121 ctags -R --languages=sh --kinds-sh=f --regex-sh='/^[ \t]*function[ \t]+([A-Za-z0-9_][A-Za-z0-9_-]*)[ \t]*\(\)?[ \t]*\{/\1/f,function/' --regex-sh='/^[ \t]*([A-Za-z0-9_][A-Za-z0-9_-]*)[ \t]*\(\)[ \t]*\{/\1/f,function/' -f .tags .
# 2122 v ~/.ctags.d/sh-hyphen-funcs.conf
ctags -R -f .tags .
