/* Original project UI. No external requests, analytics or dependencies. */
(() => {
  const root = document.documentElement;
  const lang = document.body.dataset.lang;
  const themeButton = document.querySelector('#site-theme');
  const names = lang === 'zh' ? ['主题：系统', '主题：浅色', '主题：深色'] : ['Theme: system', 'Theme: light', 'Theme: dark'];
  const states = ['system', 'light', 'dark'];
  let state = 'system';
  try { state = localStorage.getItem('ggchord-theme') || 'system'; } catch (_) { /* storage may be disabled */ }
  if (!states.includes(state)) state = 'system';
  const apply = () => {
    if (state === 'system') delete root.dataset.theme;
    else root.dataset.theme = state;
    themeButton.textContent = names[states.indexOf(state)];
    themeButton.setAttribute('aria-label', names[states.indexOf(state)]);
  };
  apply();
  themeButton.addEventListener('click', () => {
    state = states[(states.indexOf(state) + 1) % states.length];
    try { localStorage.setItem('ggchord-theme', state); } catch (_) { /* keep session theme */ }
    apply();
  });
  const languageLink = document.querySelector('.language');
  const languageBase = languageLink.href;
  const syncLanguageAnchor = () => { languageLink.href = languageBase + location.hash; };
  syncLanguageAnchor();
  addEventListener('hashchange', syncLanguageAnchor);
  const menu = document.querySelector('.navigation');
  if (matchMedia('(max-width:760px)').matches) menu.open = false;
  const input = document.querySelector('#site-search');
  const results = document.querySelector('#site-results');
  const indexURL = new URL(document.body.dataset.search, location.href);
  let index;
  let loading;
  const render = async () => {
    const query = input.value.trim().toLocaleLowerCase();
    results.replaceChildren();
    if (!query) return;
    try {
      loading ||= fetch(indexURL).then(response => {
        if (!response.ok) throw new Error('index unavailable');
        return response.json();
      });
      index ||= await loading;
      if (input.value.trim().toLocaleLowerCase() !== query) return;
      const terms = query.split(/\s+/);
      const matches = index.filter(page => terms.every(term => (page.title + ' ' + page.text).toLocaleLowerCase().includes(term)));
      matches.sort((a, b) => Number(b.title.toLocaleLowerCase().includes(query)) - Number(a.title.toLocaleLowerCase().includes(query)));
      const count = document.createElement('p');
      count.className = 'search-count';
      count.textContent = matches.length ? `${matches.length} ${lang === 'zh' ? '页匹配，最多显示 12 页' : 'pages found; showing up to 12'}` : document.body.dataset.empty;
      results.append(count);
      for (const page of matches.slice(0, 12)) {
        const a = document.createElement('a');
        a.href = new URL(page.url, indexURL).href;
        a.textContent = page.title;
        results.append(a);
      }
    } catch (_) {
      loading = undefined;
      results.textContent = lang === 'zh' ? '搜索暂不可用，请使用浏览目录。' : 'Search unavailable; use the navigation.';
    }
  };
  input.addEventListener('input', render);
  input.addEventListener('keydown', event => {
    if (event.key === 'Escape') { input.value = ''; results.replaceChildren(); }
  });
})();
