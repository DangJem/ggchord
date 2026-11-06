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
    const label=details.dataset.codeKind==='output'?(zh?'复制输出':'Copy output'):(zh?'复制本步':'Copy step');button.textContent=label;button.setAttribute('aria-live','polite');
    details.querySelector('summary').after(button);button.addEventListener('click',()=>copy(button,code.textContent,label));
  });
  document.querySelectorAll('.copy-full-example').forEach(button=>{
    const label=button.textContent;button.setAttribute('aria-live','polite');
    button.addEventListener('click',async()=>{
      try {const response=await fetch(button.dataset.source);if(!response.ok)throw Error('unavailable');await copy(button,await response.text(),label);}
      catch (_) {button.textContent=zh?'脚本暂不可用，请使用下载链接':'Script unavailable; use the download link';}
    });
  });
  // The search dialog is local, structured, keyboard-first and lazy-loaded.
  const dialog=document.querySelector('#site-search-dialog'),trigger=document.querySelector('#site-search-trigger');
  const input=document.querySelector('#site-search'),results=document.querySelector('#site-results'),status=document.querySelector('#search-status');
  if(dialog&&trigger&&input&&results) {
    const labels=zh?{all:'全部',course:'课程',function:'函数',parameter:'参数',task:'任务',faq:'排错'}:{all:'All',course:'Courses',function:'Functions',parameter:'Parameters',task:'Tasks',faq:'FAQ'};
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
      if(!query){status.textContent=zh?'输入任务、函数或参数。试试“酶切引线”或 link_branch。':'Search a task, function or parameter. Try “primer labels” or link_branch.';input.setAttribute('aria-expanded','false');return;}
      status.textContent=zh?'正在读取本地索引…':'Loading the local index…';
      try {
        loading ||= fetch(indexURL).then(r=>{if(!r.ok)throw Error('unavailable');return r.json();});records ||= await loading;
        if(input.value.trim().toLocaleLowerCase()!==query)return;
        const terms=tokens(query);
        const ranked=records.filter(r=>filter==='all'||(r.kind||'course')===filter).map(r=>{
          const title=r.title.toLocaleLowerCase(),keywords=(r.keywords||'').toLocaleLowerCase(),body=(r.text||'').toLocaleLowerCase();
          if(!terms.every(t=>(title+' '+keywords+' '+body).includes(t)))return null;
          let score=title===query?1000:title.includes(query)?150:0;
          score+=terms.reduce((n,t)=>n+(title.includes(t)?20:0)+(keywords.includes(t)?15:0),0);
          if(r.kind==='parameter'||r.kind==='task')score+=8;
          return {r,score};
        }).filter(Boolean).sort((a,b)=>b.score-a.score);
        status.textContent=ranked.length?(zh?`${ranked.length} 项匹配，显示前 20 项`:`${ranked.length} matches; showing up to 20`):(zh?'没有找到结果。试试函数全名、另一个任务词，或使用任务索引。':'No results. Try a full function name, another task, or the task index.');
        input.setAttribute('aria-expanded',String(ranked.length>0));
        for(const {r} of ranked.slice(0,20)) {
          const li=document.createElement('li');li.setAttribute('role','presentation');
          const a=document.createElement('a');a.href=new URL(r.url,indexURL).href;a.id='search-result-'+items.length;a.setAttribute('role','option');a.tabIndex=-1;a.setAttribute('aria-selected','false');
          const kind=document.createElement('span');kind.className='search-kind';kind.textContent=labels[r.kind||'course']+(r.section?' · '+r.section:'');
          const title=document.createElement('strong');highlight(title,r.title,terms);
          const excerpt=document.createElement('small');const text=(r.summary||r.text||'').replace(/\s+/g,' ');const pos=Math.max(0,text.toLocaleLowerCase().indexOf(terms.find(t=>text.toLocaleLowerCase().includes(t))||'')-45);highlight(excerpt,(pos?'…':'')+text.slice(pos,pos+180)+(text.length>pos+180?'…':''),terms);
          a.append(kind,title,excerpt);li.append(a);results.append(li);items.push(a);
        }
      } catch (_) {loading=undefined;status.textContent=zh?'搜索暂不可用，请使用任务索引。':'Search unavailable; use the task index.';}
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
  function openLinkedSection(){const target=document.getElementById(decodeURIComponent(location.hash.slice(1)));let section=target?.closest('details')||(target?.nextElementSibling?.tagName==='DETAILS'?target.nextElementSibling:null);if(section&&!section.open){section.open=true;target.scrollIntoView();}}
  openLinkedSection();addEventListener('hashchange',openLinkedSection);
})();
