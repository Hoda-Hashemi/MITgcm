const pendingCaseIds = new Set(["testcase4", "testcase7"]);
const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

const navToggle = document.querySelector("[data-nav-toggle]");
const navClose = document.querySelector("[data-nav-close]");
const navBackdrop = document.querySelector("[data-nav-backdrop]");
const navLinks = Array.from(document.querySelectorAll(".side-nav a[href^='#']"));
const sections = navLinks
  .map((link) => document.querySelector(link.getAttribute("href")))
  .filter(Boolean);

function firstText(root, selector) {
  const node = root.querySelector(selector);
  return node ? node.textContent.trim() : "";
}

function typesetMath(root) {
  if (!window.MathJax || !window.MathJax.typesetPromise) {
    return;
  }

  if (window.MathJax.typesetClear) {
    window.MathJax.typesetClear([root]);
  }
  window.MathJax.typesetPromise([root]).catch(() => {});
}

function makeChip(text, className) {
  const chip = document.createElement("span");
  chip.className = `status-chip ${className}`;
  chip.textContent = text;
  return chip;
}

function appendLinkBadge(link) {
  if (link.querySelector(".nav-status")) {
    return;
  }
  const badge = document.createElement("span");
  badge.className = "nav-status";
  badge.textContent = "Pending";
  link.append(badge);
}

function collapseDetailBlock(block) {
  if (block.matches("details")) {
    return;
  }

  const heading = block.querySelector("h3");
  const copy = block.querySelector(".description-copy");
  if (!heading || !copy) {
    return;
  }

  const details = document.createElement("details");
  details.className = `${block.className} details-panel compact-details`;
  details.dataset.compactDetail = "";

  const summary = document.createElement("summary");
  summary.textContent = heading.textContent.trim();

  details.append(summary, copy);
  block.replaceWith(details);
}

function decorateExperiment(section) {
  const topCopy = section.querySelector(".experiment-top > div:first-child");
  if (!topCopy) {
    return;
  }

  if (!topCopy.querySelector(".experiment-summary")) {
    const summaryText =
      firstText(section, ".detail-definition .description-copy p") ||
      firstText(section, ".detail-expected .description-copy p");
    if (summaryText) {
      const summary = document.createElement("p");
      summary.className = "experiment-summary";
      summary.textContent = summaryText;
      topCopy.append(summary);
    }
  }

  const summaryText = firstText(topCopy, ".experiment-summary").toLowerCase();
  const hasNoMedia = Boolean(section.querySelector(".diagnostics-stack .empty"));
  const isPending =
    pendingCaseIds.has(section.id) ||
    summaryText.includes("scaffold") ||
    (hasNoMedia && summaryText.includes("requires"));

  section.classList.toggle("is-pending", isPending);

  if (isPending && !topCopy.querySelector(".status-chip")) {
    const label = topCopy.querySelector(".section-label");
    const chip = makeChip("Scaffold / pending validation", "status-pending");
    if (label) {
      label.insertAdjacentElement("afterend", chip);
    } else {
      topCopy.prepend(chip);
    }
  }

  section.querySelectorAll(".metric-card").forEach((card) => {
    const label = firstText(card, "span").toLowerCase();
    const value = card.querySelector("strong");
    if (!value) {
      return;
    }
    if (label === "equations") {
      value.textContent = "Expandable";
    }
    if (isPending && label === "primary") {
      value.textContent = "Pending validation";
    }
  });
}

function decoratePage() {
  document.querySelectorAll(".detail-equations, .detail-parameters").forEach(collapseDetailBlock);

  document.querySelectorAll(".experiment-section").forEach((section) => {
    decorateExperiment(section);
  });

  pendingCaseIds.forEach((id) => {
    document.querySelectorAll(`a[href="#${id}"]`).forEach((link) => {
      link.classList.add("is-pending-link");
      appendLinkBadge(link);
    });
  });

  document.querySelectorAll(".diagnostics-stack .empty").forEach((empty) => {
    empty.textContent = "Scaffold / pending validation. No validated media linked yet.";
  });
}

function setNavOpen(open) {
  document.body.classList.toggle("nav-open", open);
  if (navToggle) {
    navToggle.setAttribute("aria-expanded", String(open));
  }
}

function keepActiveLinkVisible(activeLink) {
  const nav = activeLink.closest(".side-nav");
  if (!nav) {
    return;
  }

  const navRect = nav.getBoundingClientRect();
  const linkRect = activeLink.getBoundingClientRect();
  const gutter = 16;

  if (linkRect.top < navRect.top + gutter) {
    nav.scrollTop -= navRect.top + gutter - linkRect.top;
  } else if (linkRect.bottom > navRect.bottom - gutter) {
    nav.scrollTop += linkRect.bottom - navRect.bottom + gutter;
  }
}

function setActiveNav(hash) {
  let activeLink = null;
  navLinks.forEach((link) => {
    const active = link.getAttribute("href") === hash;
    link.classList.toggle("is-active", active);
    link.classList.toggle("active", active);
    if (active) {
      activeLink = link;
    }
  });

  document.querySelectorAll(".nav-group").forEach((group) => {
    const hasActive = activeLink ? group.contains(activeLink) : false;
    group.classList.toggle("has-active", hasActive);
    if (hasActive) {
      group.open = true;
    }
  });

  if (activeLink) {
    keepActiveLinkVisible(activeLink);
  }
}

function scrollToSection(hash) {
  const target = document.querySelector(hash);
  if (!target) {
    return;
  }
  target.scrollIntoView({ block: "start", behavior: reduceMotion ? "auto" : "smooth" });
  window.history.replaceState(null, "", hash);
  setActiveNav(hash);
}

decoratePage();
window.addEventListener("load", () => typesetMath(document.body));

if (navToggle) {
  navToggle.addEventListener("click", () => {
    setNavOpen(!document.body.classList.contains("nav-open"));
  });
}

[navClose, navBackdrop].forEach((control) => {
  if (control) {
    control.addEventListener("click", () => setNavOpen(false));
  }
});

navLinks.forEach((link) => {
  link.addEventListener("click", (event) => {
    const hash = link.getAttribute("href");
    if (!hash || hash === "#") {
      return;
    }
    event.preventDefault();
    scrollToSection(hash);
    setNavOpen(false);
  });
});

let activePending = false;

function updateActiveFromScroll() {
  const anchorOffset = Math.min(180, window.innerHeight * 0.28);
  let current = sections[0];
  sections.forEach((section) => {
    if (section.getBoundingClientRect().top <= anchorOffset) {
      current = section;
    }
  });
  if (current) {
    setActiveNav(`#${current.id}`);
  }
  activePending = false;
}

function requestActiveUpdate() {
  if (activePending) {
    return;
  }
  activePending = true;
  window.requestAnimationFrame(updateActiveFromScroll);
}

if (sections.length) {
  updateActiveFromScroll();
  window.addEventListener("scroll", requestActiveUpdate, { passive: true });
  window.addEventListener("resize", requestActiveUpdate);
}

function closeModal(modal) {
  const image = modal.querySelector("[data-modal-img]");
  modal.hidden = true;
  if (image) {
    image.removeAttribute("src");
  }
}

document.addEventListener("click", (event) => {
  const tab = event.target.closest("[data-tab-target]");
  if (tab) {
    const tabList = tab.closest(".tab-list");
    const tabs = tabList ? tabList.querySelectorAll("[data-tab-target]") : [];
    const tabsRoot = tab.closest(".field-tabs");
    const panels = tabsRoot ? tabsRoot.querySelectorAll("[data-tab-panel]") : [];
    const target = tab.getAttribute("data-tab-target");

    tabs.forEach((button) => {
      const active = button === tab;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-selected", String(active));
    });

    panels.forEach((panel) => {
      const active = panel.getAttribute("data-tab-panel") === target;
      panel.classList.toggle("is-active", active);
      panel.hidden = !active;
    });

    if (tabsRoot) {
      typesetMath(tabsRoot);
    }
    return;
  }

  const openButton = event.target.closest("[data-modal-src]");
  const modal = document.querySelector("[data-modal]");
  if (openButton && modal) {
    const image = modal.querySelector("[data-modal-img]");
    const caption = modal.querySelector("[data-modal-caption]");
    image.src = openButton.getAttribute("data-modal-src");
    image.alt = openButton.getAttribute("data-modal-caption") || "";
    caption.textContent = openButton.getAttribute("data-modal-caption") || "";
    modal.hidden = false;
    return;
  }

  const detailLink = event.target.closest("a[href^='#']");
  if (detailLink) {
    const target = document.querySelector(detailLink.getAttribute("href"));
    if (target instanceof HTMLDetailsElement) {
      target.open = true;
      typesetMath(target);
    }
  }

  if (modal && (event.target.closest("[data-modal-close]") || event.target === modal)) {
    closeModal(modal);
  }
});

document.addEventListener(
  "toggle",
  (event) => {
    const details = event.target;
    if (
      details instanceof HTMLDetailsElement &&
      details.open &&
      (details.matches(".compact-details") || details.querySelector(".math-line"))
    ) {
      typesetMath(details);
    }
  },
  true
);

document.addEventListener("keydown", (event) => {
  if (event.key !== "Escape") {
    return;
  }

  if (document.body.classList.contains("nav-open")) {
    setNavOpen(false);
    return;
  }

  const modal = document.querySelector("[data-modal]");
  if (modal && !modal.hidden) {
    closeModal(modal);
  }
});
