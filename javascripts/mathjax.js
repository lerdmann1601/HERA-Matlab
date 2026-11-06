window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    processEscapes: true,
    processEnvironments: true
  }
};

document$.subscribe(() => { 
  MathJax.typesetPromise()
})
