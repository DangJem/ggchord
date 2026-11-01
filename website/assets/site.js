/* Local search and paired-language navigation; Bootstrap/pkgdown owns theme/UI. */
(() => {
  const lang = document.body.dataset.lang;
  const link = document.querySelector('#site-language');
  if (link) {
    const base = link.href;
    const sync = () => { link.href = base + location.hash; };
    sync();
    addEventListener('hashchange', sync);
    link.addEventListener('click', () => {
      try { localStorage.setItem('ggchord-language', lang === 'en' ? 'zh' : 'en'); } catch (_) { /* optional */ }
    });
  }
  const input = document.querySelector('#site-search');
  const results = document.querySelector('#site-results');
  if (!input || !results) return;
  const indexURL = new URL(document.body.dataset.search, location.href);
  let index, loading;
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
      matches.sort((a,b) => Number(b.title.toLocaleLowerCase().includes(query)) - Number(a.title.toLocaleLowerCase().includes(query)));
      const count = document.createElement('p');
      count.className = 'search-count';
      count.textContent = matches.length ? `${matches.length} ${lang === 'zh' ? '页匹配，最多显示 12 页' : 'pages found; showing up to 12'}` : (lang === 'zh' ? '未找到匹配内容。' : 'No matching pages.');
      results.append(count);
      for (const page of matches.slice(0,12)) {
        const a = document.createElement('a');
        a.href = new URL(page.url,indexURL).href;
        a.textContent = page.title;
        results.append(a);
      }
    } catch (_) {
      loading = undefined;
      results.textContent = lang === 'zh' ? '搜索暂不可用，请使用顶部目录。' : 'Search unavailable; use the navigation.';
    }
  };
  input.addEventListener('input', render);
  input.addEventListener('keydown', event => {
    if (event.key === 'Escape') { input.value = ''; results.replaceChildren(); }
  });
})();
