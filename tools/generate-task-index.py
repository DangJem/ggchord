#!/usr/bin/env python3
"""Derive paired task indexes from a single reviewed task catalogue."""
import json
from documentation_lib import ROOT, WEB

def main():
    tasks=json.loads((WEB/'teaching-tasks.json').read_text())
    exports=set(__import__('re').findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),__import__('re').M))
    aliases={alias:p.stem for p in (ROOT/'man').glob('*.Rd') for alias in __import__('re').findall(r'\\alias\{([^}]+)\}',p.read_text())}
    for lang in ('en','zh'):
        han=lang=='zh'
        lines=[('按想做的事情选择课程。A 为弦图，B 为共线性图，C 为单序列圆图；带条件的条目先查看[兼容边界](compatibility.html)。基础内容完成后，再按需选择进阶内容。可用浏览器查找任务关键词，站内搜索也会索引本页。' if han else 'Choose a lesson by task. A means chord, B collinear, C single-sequence circular. Check [compatibility boundaries](compatibility.html) for conditional entries. Complete the basics, then select advanced topics as needed. Browser find and site search index these tasks.'),
          '\n## '+('功能与参数索引' if han else 'Tasks and parameters')+' {#tasks}\n',
          '| '+('任务 | 函数 | 关键参数 | 模式 | 课程与层级' if han else 'Task | Function | Key parameters | Mode | Lesson and level')+' |',
          '| --- | --- | --- | --- | --- |']
        for t in tasks:
            assert t['api'] in exports,t['api']
            slug,_,anchor=t['page'].partition('#')
            assert (WEB/lang/(slug+'.md')).exists(),slug
            target=slug+'.html'+('#'+anchor if anchor else '')
            level={'basics':('基础' if han else 'Basics'),'advanced':('进阶' if han else 'Advanced')}[t['level']]
            lines.append(f"| {t['task'][lang]} | [`{t['api']}()`](reference/{aliases[t['api']]}.html) | `{t['parameters']}` | {t['modes']} | [{level}]({target}) |")
        lines+=['\n## '+('三条独立学习路线' if han else 'Three independent learning paths')+' {#paths}\n',
         ('弦图：骨架 → 比对与基因 → 布局 → 连接路径 → 标签 → 四基因组案例。共线性图：骨架 → 输入与比对 → 行距 → 对齐与方向 → 跨行与局部窗口 → 四基因组案例。圆图：骨架 → 特征与标签 → 样式 → 注释布局 → 引物与酶切 → 轨道 → 完整 pUC19 → 可选参考图库。' if han else 'Chord: backbones → links and genes → geometry → routing → labels → four-genome case. Collinear: backbones → input and links → spacing → alignment and orientation → cross-row links and windows → four-genome case. Circular: backbone → features and labels → style → annotation layout → primers and restriction sites → tracks → complete pUC19 → optional reference gallery.'),
         '\n'+('固定样式、映射、标度、图例和主题统一在[配色与图例](styling.html)主讲；[输出](output.html)集中说明预览、尺寸和保存。函数参考提供完整签名，任务表提供教学入口，并不表示所有参数已经达到同样教学深度。' if han else '[Mappings, scales and guides](styling.html) is the shared home for appearance concepts; [output](output.html) covers preview, dimensions and saving. Reference pages provide full signatures; task links locate teaching and do not certify equal depth for every parameter.')]
        target=WEB/lang/'task-index.md'
        result='\n'.join(lines)+'\n'
        if '--check' in __import__('sys').argv:
            assert target.read_text()==result, 'Stale task index: '+lang
        else:
            target.write_text(result)
if __name__=='__main__':main()
