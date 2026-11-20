(() => {
  const button = document.getElementById('api-theme');
  if (!button) return;
  const media = window.matchMedia('(prefers-color-scheme: dark)');
  function setTheme(value) {
    document.documentElement.dataset.bsTheme = value;
    button.setAttribute('aria-pressed', value === 'dark' ? 'true' : 'false');
  }
  let saved;
  try { saved = localStorage.getItem('ggchord-api-theme'); } catch (_) { /* Optional preference. */ }
  setTheme(saved === 'dark' || saved === 'light' ? saved : media.matches ? 'dark' : 'light');
  button.hidden = false;
  button.addEventListener('click', () => {
    const theme = document.documentElement.dataset.bsTheme === 'dark' ? 'light' : 'dark';
    setTheme(theme);
    try { localStorage.setItem('ggchord-api-theme', theme); } catch (_) { /* Optional preference. */ }
  });
})();
