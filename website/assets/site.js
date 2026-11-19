/* Progressive static documentation UI. Bootstrap/pkgdown owns theme selection. */
(() => {
  'use strict';
  const zh = document.documentElement.lang.startsWith('zh');
  const language = document.querySelector('.language-control');
  const otherLanguage = document.querySelector('#site-language');
  if (otherLanguage) {
    const base=otherLanguage.href.split('#')[0];
    const sync=()=>{otherLanguage.href=base+location.hash;}; sync();addEventListener('hashchange',sync);
  }
  language?.querySelectorAll('a[data-language]').forEach(a=>a.addEventListener('click',()=>{
    try {localStorage.setItem('ggchord-language',a.dataset.language);} catch (_) { /* optional preference */ }
  }));
  document.addEventListener('click',e=>{if(language&&!language.contains(e.target))language.open=false;});
  language?.addEventListener('keydown',e=>{if(e.key==='Escape'){language.open=false;language.querySelector('summary').focus();}});
  const header=document.querySelector('.navbar');
  if(header&&'ResizeObserver'in window)new ResizeObserver(entries=>{document.documentElement.style.setProperty('--header-height',Math.ceil(entries[0].target.getBoundingClientRect().height)+'px');}).observe(header);
  // Copy actions announce exactly which range was copied.
  const copy=async(button,text,label)=>{
    let ok=false;
    try {await navigator.clipboard.writeText(text);ok=true;} catch (_) {
      const area=document.createElement('textarea');area.value=text;area.readOnly=true;
      area.style.cssText='position:fixed;left:0;top:0;opacity:0';document.body.append(area);area.select();
      try {ok=document.execCommand('copy');} catch (_) { /* manual fallback */ } area.remove();button.focus();
    }
    button.textContent=ok?(zh?'已复制':'Copied'):(zh?'复制失败，请下载脚本':'Copy failed; download the script');
    setTimeout(()=>{button.textContent=label;},2200);
  };
  document.querySelectorAll('details.code-disclosure').forEach(details=>{
    const code=details.querySelector('pre code')||details.querySelector('pre');if(!code)return;
    const button=document.createElement('button');button.type='button';button.className='copy-code';
    const kind = details.closest('[data-example-kind]')?.dataset.exampleKind;
    const label = kind === 'plot' ? (zh ? '复制绘图代码' : 'Copy plot code') :
      details.dataset.codeKind === 'signature' ? (zh ? '复制调用格式' : 'Copy usage syntax') :
      details.dataset.codeKind === 'output' ? (zh ? '复制输出' : 'Copy output') : (zh ? '复制代码' : 'Copy code');
    button.textContent=label;button.setAttribute('aria-live','polite');
    details.querySelector('summary').after(button);
    button.addEventListener('click', () => {
      const inputLines = [...code.querySelectorAll('span.r-in')];
      const text = inputLines.length ? inputLines.map(line => line.textContent).join('\n') : code.textContent;
      copy(button, text, label);
    });
  });

  // Filter exact API/parameter names without replacing the full reference.
  const parameterFilter=document.querySelector('#parameter-filter');
  if(parameterFilter){
    const entries=[...document.querySelectorAll('main details')].filter(d=>d.querySelector(':scope > summary code'));
    const notice=document.querySelector('#parameter-filter-status');
    const update=()=>{
      const query=parameterFilter.value.trim().toLocaleLowerCase();let count=0;
      entries.forEach(details=>{
        const api=details.querySelector(':scope > summary code').textContent.toLocaleLowerCase();
        const rows=[...details.querySelectorAll('tbody tr')];
        const apiMatches=api.includes(query);
        let matching=0;
        rows.forEach(row=>{const name=row.querySelector('td')?.textContent.toLocaleLowerCase()||'';row.hidden=Boolean(query)&&!apiMatches&&!name.includes(query);if(!row.hidden)matching++;});
        details.hidden=Boolean(query)&&!apiMatches&&!matching;
        if(!details.hidden){count++;if(query)details.open=true;}
      });
      notice.textContent=zh?`${count} 个函数匹配`:`${count} functions match`;
    };
    parameterFilter.addEventListener('input',update);
  }
  // The search dialog is local, structured, keyboard-first and lazy-loaded.
  const dialog=document.querySelector('#site-search-dialog'),trigger=document.querySelector('#site-search-trigger');
  const shortcut=trigger?.querySelector('kbd');
  if(shortcut)shortcut.textContent=/Mac|iPhone|iPad/.test(navigator.platform)?'⌘ K':'Ctrl K';
  const input=document.querySelector('#site-search'),results=document.querySelector('#site-results'),status=document.querySelector('#search-status');
  if(dialog&&trigger&&input&&results) {
    const labels=zh?{all:'全部',guide:'指南',function:'函数',parameter:'参数',task:'任务',tool:'查询工具',faq:'排错'}:{all:'All',guide:'Guides',function:'Functions',parameter:'Parameters',task:'Tasks',tool:'Tools',faq:'FAQ'};
    const indexURL=new URL(document.body.dataset.search,location.href);
    let records,loading,filter='all',selected=-1,items=[],invoker;
    const tokens=q=>q.match(/[a-z_][a-z0-9_.]*|[\p{Script=Han}]+/giu)?.flatMap(s=>/^[\p{Script=Han}]{3,}$/u.test(s)?(s.match(/.{1,2}/gu)||[]):[s])||[];
    const highlight=(element,value,terms)=>{
      const fragments=terms.filter(Boolean).map(t=>t.replace(/[.*+?^${}()|[\]\\]/g,'\\$&'));
      if(!fragments.length){element.textContent=value;return;}
      const re=new RegExp('('+fragments.join('|')+')','ig');
      value.split(re).forEach((part,i)=>{const node=i%2?document.createElement('mark'):document.createTextNode(part);if(i%2)node.textContent=part;element.append(node);});
    };
    const select=n=>{
      if(!items.length)return;selected=(n+items.length)%items.length;
      items.forEach((item,i)=>item.setAttribute('aria-selected',String(i===selected)));
      input.setAttribute('aria-activedescendant',items[selected].id);items[selected].scrollIntoView({block:'nearest'});
    };
    const render=async()=>{
      const query=input.value.trim().toLocaleLowerCase();results.replaceChildren();items=[];selected=-1;input.removeAttribute('aria-activedescendant');
      if(!query){status.textContent=zh?'输入任务、函数或参数，例如「覆盖度」或 geom_feature。':'Search for a task, function or parameter, such as coverage or geom_feature.';input.setAttribute('aria-expanded','false');return;}
      status.textContent=zh?'正在读取本地索引…':'Loading the local index…';
      try {
        loading ||= fetch(indexURL).then(r=>{if(!r.ok)throw Error('unavailable');return r.json();});records ||= await loading;
        if(input.value.trim().toLocaleLowerCase()!==query)return;
        const terms=tokens(query);
        const ranked=records.filter(r=>filter==='all'||(r.kind||'guide')===filter).map(r=>{
          const title=r.title.toLocaleLowerCase(),keywords=(r.keywords||'').toLocaleLowerCase(),body=(r.text||'').toLocaleLowerCase();
          if(!terms.every(t=>(title+' '+keywords+' '+body).includes(t)))return null;
          let score=title===query?1000:title.includes(query)?150:0;
          score+=terms.reduce((n,t)=>n+(title.includes(t)?20:0)+(keywords.includes(t)?15:0),0);
          if(r.kind==='guide')score+=100;
          if(r.kind==='parameter')score-=40;
          if(r.kind==='function'&&title===query)score+=1800;
          if(r.kind==='parameter'&&title.split(' · ').at(-1)===query)score+=1200;
          if((r.queries||[]).some(q=>q.toLocaleLowerCase()===query))score+=2200;
          return {r,score};
        }).filter(Boolean).sort((a,b)=>b.score-a.score);
        const seenPages=new Set();
        const grouped=ranked.filter(({r})=>{const key=r.kind==='parameter'?r.url:r.url.split('#')[0];if(seenPages.has(key))return false;seenPages.add(key);return true;});
        status.textContent=grouped.length?(zh?`${grouped.length} 项匹配，显示 ${Math.min(20,grouped.length)} 项`:`${grouped.length} matches; showing ${Math.min(20,grouped.length)}`):(zh?'没有找到结果。试试函数全名或另一个任务词。':'No results. Try a full function name or another task.');
        input.setAttribute('aria-expanded',String(grouped.length>0));
        for(const {r} of grouped.slice(0,20)) {
          const li=document.createElement('li');li.setAttribute('role','presentation');
          const a=document.createElement('a');a.href=new URL(r.url,indexURL).href;a.id='search-result-'+items.length;a.setAttribute('role','option');a.tabIndex=-1;a.setAttribute('aria-selected','false');
          const kind=document.createElement('span');kind.className='search-kind';kind.textContent=labels[r.kind||'guide']+(r.section?' · '+r.section:'');
          const title=document.createElement('strong');highlight(title,r.title,terms);
          const excerpt=document.createElement('small');const text=(r.summary||r.text||'').replace(/\s+/g,' ');const pos=Math.max(0,text.toLocaleLowerCase().indexOf(terms.find(t=>text.toLocaleLowerCase().includes(t))||'')-45);highlight(excerpt,(pos?'…':'')+text.slice(pos,pos+180)+(text.length>pos+180?'…':''),terms);
          a.append(kind,title,excerpt);li.append(a);results.append(li);items.push(a);
        }
      } catch (_) {loading=undefined;status.textContent=zh?'搜索暂不可用，请使用绘图指南。':'Search unavailable. Use the plotting guide.';}
    };
    const filters=dialog.querySelector('.search-filters');
    Object.entries(labels).forEach(([key,label])=>{const b=document.createElement('button');b.type='button';b.textContent=label;b.setAttribute('aria-pressed',String(key===filter));filters.append(b);b.addEventListener('click',()=>{filter=key;[...filters.children].forEach(x=>x.setAttribute('aria-pressed',String(x===b)));render();});});
    let previousOverflow;const open=()=>{invoker=document.activeElement;previousOverflow=document.body.style.overflow;document.body.style.overflow='hidden';dialog.showModal();input.focus();render();};
    trigger.addEventListener('click',open);dialog.querySelector('.search-close').addEventListener('click',()=>dialog.close());
    dialog.addEventListener('close',()=>{document.body.style.overflow=previousOverflow;input.setAttribute('aria-expanded','false');invoker?.focus();});
    dialog.addEventListener('click',e=>{if(e.target===dialog){const r=dialog.getBoundingClientRect();if(e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom)dialog.close();}});
    input.addEventListener('input',render);input.addEventListener('keydown',e=>{
      if(e.key==='ArrowDown'||e.key==='ArrowUp'){e.preventDefault();select(selected+(e.key==='ArrowDown'?1:-1));}
      else if(e.key==='Enter'&&selected>=0){e.preventDefault();items[selected].click();}
    });
    document.addEventListener('keydown',e=>{const editable=/INPUT|TEXTAREA|SELECT/.test(e.target.tagName)||e.target.isContentEditable;if(((e.metaKey||e.ctrlKey)&&e.key.toLowerCase()==='k')||(!editable&&e.key==='/')){e.preventDefault();if(!document.querySelector('dialog[open]'))open();}});
    trigger.setAttribute('aria-label',zh?'搜索任务、函数和参数，快捷键 Control 或 Command K':'Search tasks, functions and parameters, Control or Command K');
  }
  const galleryControls = document.querySelector('.gallery-controls');
  if (galleryControls) {
    const groups = [...document.querySelectorAll('.gallery-family')];
    const family = galleryControls.querySelector('[name="family"]');
    const scope = galleryControls.querySelector('[name="scope"]');
    const status = galleryControls.querySelector('[role="status"]');
    function report() {
      const shown = groups.filter(group => !group.hidden).reduce((total, group) =>
        total + [...group.querySelectorAll('.gallery-card')].filter(card =>
          !card.hidden && !card.closest('details:not([open])')).length, 0);
      status.textContent = shown ? (zh ? `显示 ${shown} 张图` : `Showing ${shown} figures`) :
        (zh ? '此分类没有参数对照图。' : 'No parameter comparisons in this category.');
    }
    function apply() {
      groups.forEach(group => {
        const cards = [...group.querySelectorAll('.gallery-card')];
        cards.forEach(card => {card.hidden = scope.value === 'comparison' && card.dataset.comparison !== 'true';});
        group.hidden = (family.value !== 'all' && group.dataset.family !== family.value) ||
          cards.every(card => card.hidden);
        const more = group.querySelector('details');
        if (more) more.open = scope.value !== 'featured';
      });
      report();
    }
    function revealHash() {
      let id;
      try {id = decodeURIComponent(location.hash.slice(1));} catch (_) {return;}
      const target = document.getElementById(id);
      if (!target?.closest('.gallery-family')) return;
      family.value = 'all';
      scope.value = target.closest('.gallery-more') ? 'all' : 'featured';
      apply();
      target.scrollIntoView({block: 'start'});
    }
    groups.forEach(group => group.querySelector('details')?.addEventListener('toggle', report));
    galleryControls.hidden = false;
    galleryControls.addEventListener('change', apply);
    apply();
    revealHash();
    addEventListener('hashchange', revealHash);
  }
  // Native dialog isolates the background; full-resolution PNG pixels are unchanged.
  const figures=[...document.querySelectorAll('a.figure-link')];
  if(figures.length) {
    const viewer=document.createElement('dialog');viewer.className='figure-viewer';viewer.setAttribute('aria-labelledby','figure-viewer-title');
    viewer.innerHTML='<div class="figure-toolbar"><strong id="figure-viewer-title"></strong><button class="viewer-button figure-size" type="button"></button><a class="figure-original" target="_blank" rel="noopener"></a><a class="figure-source"></a><a class="figure-code"></a><button class="viewer-button figure-close" type="button"></button></div><div class="figure-canvas" tabindex="0"><img alt=""></div>';
    document.body.append(viewer);viewer.querySelector('.figure-canvas').setAttribute('role','region');viewer.querySelector('.figure-canvas').setAttribute('aria-label',zh?'科研图；原始大小下可用方向键、滚动或拖动查看':'Scientific figure; use arrow keys, scroll or drag at actual size');const canvas=viewer.querySelector('.figure-canvas'),image=canvas.querySelector('img'),size=viewer.querySelector('.figure-size'),close=viewer.querySelector('.figure-close'),original=viewer.querySelector('.figure-original'),source=viewer.querySelector('.figure-source'),codeLink=viewer.querySelector('.figure-code');let invoker,drag,previousOverflow;
    codeLink.textContent=zh?'对应代码':'View code';close.textContent=zh?'关闭':'Close';original.textContent=zh?'打开原始 PNG':'Open original PNG';source.textContent=zh?'来源与许可':'Sources and terms';source.href=new URL('../assets/ASSET-SOURCES.md',new URL(document.body.dataset.search,location.href)).href;
    const fit=()=>{viewer.classList.remove('native-size');size.textContent=zh?'原始大小':'Actual size';canvas.scrollTo(0,0);};
    size.addEventListener('click',()=>{const native=viewer.classList.toggle('native-size');size.textContent=native?(zh?'适应窗口':'Fit window'):(zh?'原始大小':'Actual size');});
    codeLink.addEventListener('click',()=>viewer.close());close.addEventListener('click',()=>viewer.close());viewer.addEventListener('close',()=>{document.body.style.overflow=previousOverflow;invoker?.focus();});
    canvas.addEventListener('pointerdown',e=>{if(!viewer.classList.contains('native-size')||e.pointerType==='touch')return;drag={x:e.clientX,y:e.clientY,left:canvas.scrollLeft,top:canvas.scrollTop};canvas.setPointerCapture(e.pointerId);});
    canvas.addEventListener('pointermove',e=>{if(drag){canvas.scrollLeft=drag.left-e.clientX+drag.x;canvas.scrollTop=drag.top-e.clientY+drag.y;}});canvas.addEventListener('pointerup',()=>{drag=null;});
    figures.forEach(a=>a.addEventListener('click',e=>{if(e.metaKey||e.ctrlKey||e.shiftKey||e.altKey)return;e.preventDefault();invoker=a;image.src=a.href;image.alt=a.querySelector('img')?.alt||'';viewer.querySelector('strong').textContent=image.alt;original.href=a.href;codeLink.href=a.dataset.codeUrl||location.href;fit();previousOverflow=document.body.style.overflow;document.body.style.overflow='hidden';viewer.showModal();close.focus();}));
  }
  // Right outline tracks location without moving keyboard focus.
  const links=[...document.querySelectorAll('.page-nav a[href^="#"]')],targets=links.map(a=>document.getElementById(decodeURIComponent(a.hash.slice(1)))).filter(Boolean);
  if(targets.length&&'IntersectionObserver'in window) {
    const observer=new IntersectionObserver(()=>{
      const target=targets.filter(h=>h.getBoundingClientRect().top<160).at(-1)||targets[0];
      links.forEach(a=>{if(a.hash==='#'+target.id)a.setAttribute('aria-current','location');else a.removeAttribute('aria-current');});
    },{rootMargin:'-90px 0px -55% 0px',threshold:[0,1]});targets.forEach(h=>observer.observe(h));
  }
  function openLinkedSection(){
    const target=document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if(parameterFilter&&target&&/^(param-|api-)/.test(target.id)){
      parameterFilter.value='';parameterFilter.dispatchEvent(new Event('input'));
    }
    const following=target?.nextElementSibling||target?.closest('p')?.nextElementSibling;
    const section=target?.closest('details')||(following?.tagName==='DETAILS'?following:null);
    if(section&&!section.open){section.open=true;target.scrollIntoView();}
  }
  openLinkedSection();addEventListener('hashchange',openLinkedSection);
})();

// Filter the reference index locally; the unfiltered lists work without JavaScript.
(() => {
  const input = document.getElementById('reference-filter');
  if (!input) return;
  input.closest('.reference-filter').hidden = false;
  const groups = [...document.querySelectorAll('.reference-group')];
  const count = document.getElementById('reference-filter-count');
  function filter() {
    const query = input.value.trim().toLocaleLowerCase();
    let total = 0;
    for (const group of groups) {
      const heading = group.querySelector('h2').textContent.toLocaleLowerCase();
      let matches = 0;
      for (const term of group.querySelectorAll('dt')) {
        const definition = term.nextElementSibling;
        const text = (term.textContent + ' ' + definition.textContent).toLocaleLowerCase();
        const visible = heading.includes(query) || text.includes(query);
        term.hidden = definition.hidden = !visible;
        if (visible) matches++;
      }
      group.hidden = !matches;
      total += matches;
    }
    count.textContent = !query ? '' : document.documentElement.lang.startsWith('zh')
      ? '匹配 ' + total + ' 个主题' : total + (total === 1 ? ' matching topic' : ' matching topics');
  }
  input.addEventListener('input', filter);
  window.addEventListener('hashchange', () => {
    if (groups.some(group => group.hidden && location.hash.slice(1) === group.querySelector('h2').id)) {
      input.value = '';
      filter();
      document.getElementById(location.hash.slice(1)).scrollIntoView();
    }
  });
})();
