/* Local search and paired-language navigation; Bootstrap/pkgdown owns theme/UI. */
(() => {
  const lang = document.body.dataset.lang;
  const language = document.querySelector('.language-control');
  const link = document.querySelector('#site-language');
  if (link) {
    const base = link.href.split('#')[0];
    const sync = () => { link.href = base + location.hash; };
    sync();
    addEventListener('hashchange', sync);
  }
  if (language) {
    for (const option of language.querySelectorAll('a[data-language]')) {
      option.addEventListener('click', () => {
        try { localStorage.setItem('ggchord-language', option.dataset.language); } catch (_) { /* optional */ }
      });
    }
    document.addEventListener('click', event => {
      if (!language.contains(event.target)) language.open = false;
    });
    language.addEventListener('keydown', event => {
      if (event.key === 'Escape') { language.open = false; language.querySelector('summary').focus(); }
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

/* Progressive image viewer; plain links remain usable without JavaScript. */
(() => {
  const zh = document.documentElement.lang.startsWith('zh');
  const links = [...document.querySelectorAll('a.figure-link')];
  if (!links.length) return;
  const backdrop = document.createElement('div');
  backdrop.className = 'figure-backdrop';
  backdrop.hidden = true;
  backdrop.innerHTML = '<div class="figure-viewer" role="dialog" aria-modal="true" aria-labelledby="figure-viewer-title"><div class="figure-toolbar"><strong id="figure-viewer-title"></strong><button type="button" class="figure-size"></button><a class="figure-original" target="_blank" rel="noopener"></a><button type="button" class="figure-close"></button></div><div class="figure-canvas" tabindex="0"><img alt=""></div></div>';
  document.body.append(backdrop);
  const viewer = backdrop.querySelector('.figure-viewer');
  const image = viewer.querySelector('img');
  const size = viewer.querySelector('.figure-size');
  const original = viewer.querySelector('.figure-original');
  const close = viewer.querySelector('.figure-close');
  let invoker, overflow, background = [];
  original.textContent = zh ? '打开原图' : 'Open original';
  close.textContent = zh ? '关闭' : 'Close';
  const fit = () => { viewer.classList.remove('native-size'); size.textContent = zh ? '原始大小' : 'Actual size'; };
  const hide = () => {
    backdrop.hidden = true;
    document.body.style.overflow = overflow;
    for (const [element, previous] of background) {
      if (previous === null) element.removeAttribute('aria-hidden');
      else element.setAttribute('aria-hidden',previous);
    }
    background = [];
    invoker?.focus();
  };
  size.addEventListener('click', () => {
    const native = viewer.classList.toggle('native-size');
    size.textContent = native ? (zh ? '适应窗口' : 'Fit window') : (zh ? '原始大小' : 'Actual size');
  });
  close.addEventListener('click',hide);
  backdrop.addEventListener('click',event => { if (event.target === backdrop) hide(); });
  viewer.addEventListener('keydown',event => {
    if (event.key === 'Escape') { event.preventDefault(); hide(); }
    if (event.key === 'Tab') {
      const items = [...viewer.querySelectorAll('button,a[href],[tabindex="0"]')];
      const first = items[0], last = items[items.length-1];
      if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
      else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
    }
  });
  for (const link of links) link.addEventListener('click',event => {
    if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    invoker = link;
    image.src = link.href;
    image.alt = link.querySelector('img').alt;
    viewer.querySelector('#figure-viewer-title').textContent = image.alt;
    original.href = link.href;
    fit();
    background = [...document.body.children].filter(e=>e!==backdrop).map(e=>[e,e.getAttribute('aria-hidden')]);
    for (const [element] of background) element.setAttribute('aria-hidden','true');
    overflow = document.body.style.overflow;
    document.body.style.overflow = 'hidden';
    backdrop.hidden = false;
    close.focus();
  });
})();

/* Copy only executable text; disclosure labels are outside the copied range. */
(() => {
  const zh = document.documentElement.lang.startsWith('zh');
  for (const disclosure of document.querySelectorAll('details.code-disclosure')) {
    const code = disclosure.querySelector('pre code') || disclosure.querySelector('pre');
    if (!code) continue;
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'copy-code';
    const label = zh ? '复制代码' : 'Copy code';
    button.textContent = label;
    button.setAttribute('aria-label', label);
    button.setAttribute('aria-live', 'polite');
    disclosure.querySelector('summary').after(button);
    button.addEventListener('click', async () => {
      const text = code.textContent;
      let copied = false;
      try {
        if (navigator.clipboard?.writeText) {
          await navigator.clipboard.writeText(text);
          copied = true;
        }
      } catch (_) { /* Fall back to the document copy action. */ }
      if (!copied) {
        const field = document.createElement('textarea');
        field.value = text;
        field.readOnly = true;
        field.setAttribute('aria-hidden', 'true');
        field.style.cssText = 'position:fixed;opacity:0;left:0;top:0';
        document.body.append(field);
        field.select();
        try { copied = document.execCommand('copy'); } catch (_) { /* Select code below. */ }
        field.remove();
        button.focus();
      }
      if (!copied) {
        const range = document.createRange();
        range.selectNodeContents(code);
        const selection = getSelection();
        selection.removeAllRanges();
        selection.addRange(range);
      }
      button.textContent = copied ? (zh ? '已复制' : 'Copied') : (zh ? '已选中，请复制' : 'Selected; copy manually');
      setTimeout(() => { button.textContent = label; }, 2500);
    });
  }
})();

/* The left outline follows the current section; the right guide follows the task. */
(() => {
  const links = [...document.querySelectorAll('.page-nav a[href^="#"]')];
  const targets = links.map(link => document.getElementById(decodeURIComponent(link.hash.slice(1)))).filter(Boolean);
  if (!targets.length || !('IntersectionObserver' in window)) return;
  const observer = new IntersectionObserver(entries => {
    const visible = entries.filter(entry => entry.isIntersecting).sort((a,b) => a.boundingClientRect.top - b.boundingClientRect.top)[0];
    if (!visible) return;
    for (const link of links) {
      if (link.hash === '#' + visible.target.id) link.setAttribute('aria-current','location');
      else link.removeAttribute('aria-current');
    }
  }, {rootMargin: '-100px 0px -60% 0px'});
  targets.forEach(target => observer.observe(target));
})();
